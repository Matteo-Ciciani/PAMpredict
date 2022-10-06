import os
import sys
import logging
import numpy as np
import pandas as pd
import shutil
from Bio import Seq,SeqIO

class Blastn:
    
    def __init__(self, obj):
        for key, val in vars(obj).items():
            setattr(self, key, val)
            
    def unique_spacers(self):
        # save unique spacer sequences
        unique_spacer_dict = {self.spacer_dict[ID]:ID for ID in self.spacer_dict}
        self.unique_spacer_dict = {unique_spacer_dict[seq]:seq for seq in unique_spacer_dict}
        unique_spacer_records = [SeqIO.SeqRecord(Seq.Seq(self.unique_spacer_dict[ID]), id=ID,
            name='', description='') for ID in self.unique_spacer_dict]
        self.unique_spacers_file = os.path.join(self.outdir, 'unique_spacers.fna')
        logging.info('{} unique spacers read.'.format(len(unique_spacer_records)))
        if len(unique_spacer_records)<10:
            logging.warning('Less than 10 unique spacers present, predictions might be unreliable!')
        SeqIO.write(unique_spacer_records, self.unique_spacers_file, 'fasta')
    
    def align_spacers(self):
        self.blastout_files = {}
        self.blastdir = os.path.join(self.outdir, 'blastn' if self.keep_tmp else '.blastn')
        if not os.path.isdir(self.blastdir):
            os.makedirs(self.blastdir)
        for phage_fasta in self.db_fasta:
            logging.info('Running blastn on ' + phage_fasta)
            self.blastout_files[phage_fasta] = os.path.join(self.blastdir, 'blastout_{}.tsv'.format(
                ''.join(phage_fasta.split('.')[:-1])))
            blastn_command = "blastn -out {} -outfmt \"6  qseqid sseqid pident nident qlen slen length mismatch gapopen qseq sseq qstart qend sstart send evalue sstrand bitscore\" -query {} -db {} -task blastn-short -gapopen 1 -gapextend 2 -penalty -1 -reward 1 -evalue 1 -word_size 10 -num_threads {}"
            blastn_command = blastn_command.format(self.blastout_files[phage_fasta], 
                self.unique_spacers_file, os.path.join(self.phage_db, phage_fasta), self.threads)
            os.system(blastn_command) # should work
    
    def read_result(self):
        logging.info('Reading matches...')
        cols = 'qseqid sseqid pident nident qlen slen length mismatch gapopen qseq sseq qstart qend sstart send evalue sstrand bitscore'.split(' ')
        self.blast_results = {}
   
        for phage_fasta in self.db_fasta:
            self.blast_results[phage_fasta] = pd.read_csv(self.blastout_files[
                phage_fasta], sep='\t', header=None, names=cols)
    
    def filter_spacers(self):
        logging.info('Filtering matches...')
        self.filtered_matches = {}
        self.filtered_match_count = {}

        for phage_fasta in self.db_fasta:
            near_perfect_matches = self.blast_results[phage_fasta][self.blast_results[
                phage_fasta]['qlen']-self.blast_results[phage_fasta][
                'nident'] <= self.max_diff]
            full_len_matches = near_perfect_matches[near_perfect_matches[
                'length']==near_perfect_matches['qlen']]
            self.filtered_matches[phage_fasta] = full_len_matches
            for ID in full_len_matches['qseqid'].unique():
                if ID not in self.filtered_match_count:
                    self.filtered_match_count[ID] = 0
                self.filtered_match_count[ID] += full_len_matches[
                    full_len_matches['qseqid']==ID].shape[0]
    
    def write_stats(self):
        
        with open(os.path.join(self.outdir, 'spacer_alignment_stats.tsv'), 'w') as fh:
            for phage_fasta in self.db_fasta:
                fh.write('Total # of matches in {}:\t{}\n'.format(
                    phage_fasta, self.blast_results[phage_fasta].shape[0]))
            fh.write('\n')
            for phage_fasta in self.db_fasta:
                fh.write('Total # of matches after filtering in {}:\t{}\n'.format(
                    phage_fasta, self.filtered_matches[phage_fasta].shape[0]))
            fh.write('Total # of unique spacers: {}\n'.format(len(
                self.unique_spacer_dict)))
            fh.write('Total # of mapped spacers: {}\n'.format(len(
                self.filtered_match_count)))
            fh.write('\n')
            fh.write('Stats by spacer:\n\n')
            fh.write('Spacer ID\t' + '\t'.join(['# matches in '+phage_fasta+
                '\t # matches after filtering in '+phage_fasta for phage_fasta
                in self.db_fasta]) + '\n')
                
            for ID in self.unique_spacer_dict:
                fh.write(ID + ''.join(['\t{}\t{}'.format(self.blast_results[
                    phage_fasta][self.blast_results[phage_fasta]['qseqid']==ID].shape[0],
                    self.filtered_matches[phage_fasta][self.filtered_matches[
                    phage_fasta]['qseqid']==ID].shape[0]) for phage_fasta
                    in self.db_fasta]) + '\n')
        
        if len(self.filtered_match_count) < 1:
            logging.info('No full-lenght, near-perfect macth found!')
            logging.info('Consider increasing the maximum number of differences between spacers and matches with --max_diff (this might result in inaccurate predictions).')
        if len(self.filtered_match_count) < 10:
            logging.warning('Only {} mapped spacers, predictions might be inaccurate.'.format(
                len(self.filtered_match_count)))
            
    def blast_cleanup(self):
        if not self.keep_tmp:
            shutil.rmtree(self.blastdir)
        if len(self.filtered_match_count) < 1:
            sys.exit(0)
    
    def run_blastn(self):
        self.unique_spacers()
        self.align_spacers()
        self.read_result()
        self.filter_spacers()
        self.write_stats()
        self.blast_cleanup()
