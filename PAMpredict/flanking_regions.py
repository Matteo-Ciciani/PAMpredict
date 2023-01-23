import os
import sys
import logging
import numpy as np
import pandas as pd
import subprocess
import multiprocessing as mp
from pysam import FastaFile
from Bio import Seq, SeqIO
from Bio import pairwise2
from collections import Counter
import shutil

class Flanking_Regions:
    
    def __init__(self, obj):
        for key, val in vars(obj).items():
            setattr(self, key, val)
        self.flank_len = 30
            
    def read_fasta_db(self):
        logging.info('Reading flanking sequences from phage databases...')
        self.indexed_fasta = {}
        for phage_fasta in self.db_fasta:
            if phage_fasta + '.fai' not in self.all_db_files:
                logging.info('No index file found for '+phage_fasta)
                logging.info('Building index, this is run only once and might take some time...')
                subprocess.run(['samtools', 'faidx', os.path.join(self.phage_db, phage_fasta)])
                if not os.path.isfile(os.path.join(self.phage_db, phage_fasta)):
                    logging.error('Samtools faild to build the index, exiting...')
                    sys.exit(1)
            # FastaFile needs a .fai index
            self.indexed_fasta[phage_fasta] = FastaFile(os.path.join(
                self.phage_db, phage_fasta))
            
    def fetch(self):
        for phage_fasta in self.db_fasta:
            upstream_seqs, downstream_seqs = [], []
            for i in self.filtered_matches[phage_fasta].index:
                sseqid = self.filtered_matches[phage_fasta].loc[i, 'sseqid']
                sstart = self.filtered_matches[phage_fasta].loc[i, 'sstart']
                send = self.filtered_matches[phage_fasta].loc[i, 'send']
                sstrand = self.filtered_matches[phage_fasta].loc[i, 'sstrand']
                seq = self.indexed_fasta[phage_fasta].fetch(reference=sseqid)
                if sstrand=='plus':
                    upstream_seqs.append(seq[max(0,sstart-self.flank_len-1):max(0,sstart-1)])
                    downstream_seqs.append(seq[send:send+self.flank_len])
                else:
                    upstream_seqs.append(str(Seq.Seq(seq[sstart:sstart+self.flank_len]
                        ).reverse_complement()))
                    downstream_seqs.append(str(Seq.Seq(seq[max(0,
                        send-self.flank_len-1):max(0,send-1)]).reverse_complement()))
            self.filtered_matches[phage_fasta] = self.filtered_matches[
                phage_fasta].assign(**{'Upstream':upstream_seqs,
                'Downstream':downstream_seqs})
            if self.keep_tmp:
                for phage_fasta in self.filtered_matches:
                    self.filtered_matches[phage_fasta].to_csv(os.path.join(
                        self.blastdir, '.'.join(phage_fasta.split(
                        '.')[:-1])+'_filtered_matches_with_flanking_sequences.tsv'),
                        sep='\t')
        self.filtered_matches = pd.concat([self.filtered_matches[
            phage_fasta] for phage_fasta in self.filtered_matches])
        self.filtered_matches = self.filtered_matches.reset_index()
                
    def realign_one_spacer(self, spacer):
        seqs_to_align = []
        us_matches = self.filtered_matches[self.filtered_matches['qseqid']==spacer]
        # if there's only 1 sequence, just return its flanking sequences
        if us_matches.shape[0]==1:
            return us_matches.iloc[0]['Upstream'], us_matches.iloc[0]['Downstream']
        # make upstream + match + downstream sequeneces
        for i in us_matches.index:
            s = ''
            for x in [us_matches.loc[i, 'Upstream'], us_matches.loc[i,
                'sseq'], us_matches.loc[i, 'Downstream']]:
                if type(x)==str: s+=x
            seqs_to_align.append(SeqIO.SeqRecord(Seq.Seq(s), id='{}_{}'.format(spacer, i)))
        # write sequences
        fasta_to_align = os.path.join(self.realigndir, '{}_to_realign.fasta'.format(spacer))
        SeqIO.write(seqs_to_align, fasta_to_align, 'fasta')
        #align them
        aligned_fasta = os.path.join(self.realigndir, '{}_aligned.fasta'.format(spacer))
        if len(seqs_to_align) < 1000:
            os.system('muscle -quiet -align {} -output {}'.format(
                fasta_to_align, aligned_fasta))
        elif len(seqs_to_align) < 10000:
            os.system('muscle -super5 {} -output {}'.format(
                fasta_to_align, aligned_fasta))
        else:
            os.system('mafft --retree 1 --quiet --thread {} {} > {}'.format(
                self.threads, fasta_to_align, aligned_fasta))
        # read alignment
        aln = [str(rec.seq).upper() for rec in SeqIO.parse(aligned_fasta, 'fasta')]
        # compute consensus from alignment
        counters = [Counter() for i in range(len(aln[0]))]
        for i in range(len(counters)):
            for s in aln:
                counters[i][s[i]] +=1
        consensus_flank = ''.join([pd.Series(c).idxmax() for c
            in counters if len(c)>0]).replace('-', '')
        # realign spacer to consensus, defining consensus up and down
        alignment = pairwise2.align.localxs(consensus_flank,
            self.unique_spacer_dict[spacer], -1, -1)[0]
        return consensus_flank[0:alignment.start], consensus_flank[alignment.end:]
    
    def realign(self):
        logging.info('Realigning flanking regions...')
        # remove self.indexed_fasta since it's no longer useful and it messes up
        # with the parallelization (it causes a cython error as it can't be pickeld)
        self.indexed_fasta = None
        realigned_seqs = {}
        self.realigndir = os.path.join(self.outdir,
        'realignment' if self.keep_tmp else '.realignemnt')
        if not os.path.isdir(self.realigndir):
            os.makedirs(self.realigndir)
        
        # realign spacers with > 10k matches with all threads
        for ID in self.filtered_match_count:
            if self.filtered_match_count[ID]>10000:
                realigned_seqs[ID] = self.realign_one_spacer(ID)
        # realign the rest
        if self.threads > 1:
            pool = mp.Pool(self.threads)
            tmp_res = {ID:pool.apply_async(self.realign_one_spacer, [ID])
                for ID in self.filtered_match_count if ID not in realigned_seqs}
            pool.close()
            pool.join()
            for ID in tmp_res:
                realigned_seqs[ID] = tmp_res[ID].get()
        else:
            for ID in self.filtered_match_count:
                if ID not in realigned_seqs:
                    realigned_seqs[ID] = self.realign_one_spacer(ID)
        
        logging.info('Flanking regions realignment completed!')
        # reformat as dataframe
        self.consensus_flank_df = pd.DataFrame(realigned_seqs,
            index=['Upstream', 'Downstream']).transpose()
        
    def realignment_clenaup(self):
        if self.keep_tmp:
            self.consensus_flank_df.to_csv(os.path.join(self.outdir,
                'consensus_flanking_sequences.tsv'), sep='\t')
        else:
            shutil.rmtree(self.realigndir)
    
    def run_realignemnt(self):
        self.read_fasta_db()
        self.fetch()
        self.realign()
        self.realignment_clenaup()
