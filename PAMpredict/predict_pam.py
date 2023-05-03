import os
import sys
import logging
import numpy as np
from Bio import Seq,SeqIO
import pandas as pd
import numpy as np
from scipy import stats
from collections import Counter
import logomaker
import math
import matplotlib.pyplot as plt
from matplotlib import ticker

class Predict_PAM:
    
    def __init__(self, obj):
        for key, val in vars(obj).items():
            setattr(self, key, val)
            
    def make_info_df(self):
        logging.info('Generating PAM prediction...')
        self.info_df = {}
        # downstream
        self.info_df['Downstream'] = [Counter() for i in range(self.pam_len)]
        for i in range(len(self.info_df['Downstream'])):
            for seq in self.consensus_flank_df['Downstream'].values:
                if i < len(seq):
                    if seq[i]!='N':
                        self.info_df['Downstream'][i][seq[i]] +=1
        # upstream
        self.info_df['Upstream'] = [Counter() for i in range(self.pam_len)]
        for i in range(1,len(self.info_df['Upstream'])+1):
            for seq in self.consensus_flank_df['Upstream'].values:
                if i <= len(seq):
                    if seq[-i]!='N':
                        self.info_df['Upstream'][-i][seq[-i]] +=1
        # fill missing values
        for pos in self.info_df:
            for i in range(len(self.info_df[pos])):
                for base in list('ACGT'):
                    if base not in self.info_df[pos][i]:
                        self.info_df[pos][i][base] = 0
        try:
            for pos in self.info_df:
                freqs = pd.DataFrame(self.info_df[pos]).fillna(0)
                freqs = freqs.divide(freqs.sum(axis=1), axis=0).fillna(1.0/freqs.shape[1])
                if not all(freqs.sum(axis=1).apply(lambda x: math.isclose(x, 1.0))):
                    raise ValueError
                self.info_df[pos] = logomaker.transform_matrix(freqs,
                from_type='probability', to_type='information')
        except ValueError:
            logging.error('PAM frequency matrix doesn\'t sum to 1.')
            sys.exit(1)
            
    def compute_threshold(self):
        self.all_info = {pos:self.info_df[pos].sum(axis=1).values for pos in self.info_df}
        self.all_info_all_pos = [x for pos in self.all_info for x in self.all_info[pos]]
        self.threshold = min(2, max(1, pd.Series(self.all_info_all_pos).quantile(
            0.75)+1.5*stats.iqr(self.all_info_all_pos)))
            
    def compute_prediction(self):
        # check positions with more information than threshold
        self.downstream_pam = self.all_info['Downstream'][0:self.pam_len]
        self.upstream_pam = self.all_info['Upstream'][-self.pam_len:]
        self.has_downstream_pam = any(self.downstream_pam>self.threshold)
        self.has_upstream_pam = any(self.upstream_pam>self.threshold)
        
        if self.has_downstream_pam and self.has_upstream_pam:
            logging.warning('Conserved bases detected in both flanking sequences! The PAM prediction is ambiguous.')
        elif not self.has_downstream_pam and not self.has_upstream_pam:
            logging.warning('No conserved base detected in either flanking sequences!')
        else:
            logging.info('Conserved bases detected in the {} flanking sequence!'.format(self.pam_position))
            
    def reverse_pam(self, pam):
        newpam = pam[['A', 'T', 'C', 'G']][::-1]
        newpam.columns = ['T', 'A', 'G', 'C']
        newpam.index = newpam.index[::-1]
        return newpam
        
    def logo_to_seq(self):
        IUPAC = {
                ('A', 'G'): 'R',
                ('C', 'T'): 'Y',
                ('C', 'G'): 'S',
                ('A', 'T'): 'W',
                ('G', 'T'): 'K',
                ('A', 'C'): 'M',
                ('C', 'G', 'T'): 'B',
                ('A', 'G', 'T'): 'D',
                ('A', 'C', 'T'): 'H',
                ('A', 'C', 'G'): 'V'
                }
        
        if self.pam_position == 'downstream':
            info = self.info_df['Downstream']
        else:
            info = self.info_df['Upstream']
        
        # turn info df into list of letters, e.g [['N'], ['G'], ['A','G']], which is NGR
        PAM = []
        for i in info.index:
            PAM.append([])
            for base in info.columns:
                if info.loc[i, base] >= self.threshold*0.5:
                    PAM[i].append(base)
            if len(PAM[i])==0:
                PAM[i].append('N')

        # turn list of letter into string
        flat_pam = ''
        for x in PAM:
            if len(x)==1:
                flat_pam = flat_pam + x[0]
            else:
                letter = IUPAC[tuple(np.sort(x))]
                flat_pam = flat_pam + letter
            
        self.pam_seq = flat_pam.rstrip('N') if self.pam_position == 'downstream' else flat_pam.lstrip('N')
        
    def save_report(self):
        with open(os.path.join(self.outdir, 'PAM_prediction.txt'), 'w') as fh:
            fh.write('Conserved bases detected in the {} flanking sequence.\n'.format(
                self.pam_position))
            fh.write('Predicted PAM: {} (rough approxiamtion, visual inspection fo the sequence logo is recommended).\n'.format(
                self.pam_seq))
            fh.write('Inferred input spacers orientation: {}.\n'.format(
                'Reverse' if self.pam_on_wrong_strand else 'Forward'))

    def save_prediction(self):
        if self.has_downstream_pam ^ self.has_upstream_pam:
            # check if the PAM is on the wrong strand, if it is flip it 
            # (this happens when the input spacers are in the wrong orientation)
            self.pam_on_wrong_strand = (self.pam_position=='downstream' and self.has_upstream_pam) or (
                self.pam_position=='upstream' and self.has_downstream_pam)
            if self.pam_on_wrong_strand:
                self.info_df['Upstream'], self.info_df['Downstream'] = self.reverse_pam(
                    self.info_df['Downstream']), self.reverse_pam(self.info_df['Upstream'])
            # make PAM sequence
            self.logo_to_seq()
            self.save_report()
            if not self.no_plot:
                logging.info('Saving {} flanking sequence logo (PAM) and {} sequence logo (as a reference)...'.format(
                    self.pam_position, 'upstreamdownstream'.replace(self.pam_position, '')))
        else:
            if not self.no_plot:
                logging.info('Saving both flanking sequence logos...')
        
        # save PAM dataframes
        pd.DataFrame(self.info_df['Upstream'][['A','C','G','T']]).to_csv(os.path.join(self.outdir, 'upstream_flanking_sequence_info.tsv'), sep='\t')
        pd.DataFrame(self.info_df['Downstream'][['A','C','G','T']]).to_csv(os.path.join(self.outdir, 'downstream_flanking_sequence_info.tsv'), sep='\t')

    
    def plot(self):
        if not self.no_plot:
            
            # downstream
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,1.5))
            # create Logo object
            seq_logo = logomaker.Logo(self.info_df['Downstream'], ax=ax)
            # style using Logo methods
            seq_logo.style_spines(visible=False)
            seq_logo.style_spines(spines=['left', 'bottom'], visible=True)
            seq_logo.style_xticks(fmt='%d', anchor=0)
            ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
            seq_logo.ax.set_xticklabels(np.array(range(1, self.pam_len+1)), size=12)
            for y in seq_logo.ax.get_yticklabels():
                y.set_size(12)
            seq_logo.ax.set_ylabel("bits", size=15)
            ax.set_ylim(0,2)
            plt.tight_layout()
            plt.savefig(os.path.join(self.outdir, 'Downstream_flanking_sequence.'+self.format), dpi=300)
            
            # upstream
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,1.5))
            # create Logo object
            seq_logo = logomaker.Logo(self.info_df['Upstream'], ax=ax)
            # style using Logo methods
            seq_logo.style_spines(visible=False)
            seq_logo.style_spines(spines=['left', 'bottom'], visible=True)
            seq_logo.style_xticks(fmt='%d', anchor=0)
            ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
            seq_logo.ax.set_xticklabels([-x for x in range(1, self.pam_len+1)][::-1])
            for y in seq_logo.ax.get_yticklabels():
                y.set_size(12)
            seq_logo.ax.set_ylabel("bits", size=15)
            ax.set_ylim(0,2)
            plt.tight_layout()
            plt.savefig(os.path.join(self.outdir, 'Upstream_flanking_sequence.'+self.format), dpi=300)
            
    def predict(self):
        self.make_info_df()
        self.compute_threshold()
        self.compute_prediction()
        self.save_prediction()
        self.plot()
