#!/usr/bin/env python3

###############################
#       WRAPPER SCRIPT        #
###############################

import os
import sys
import argparse
import logging
from Bio import Seq,SeqIO
import pandas as pd
import numpy as np
from scipy import stats
import subprocess
import multiprocessing as mp
import logomaker
import math

from controller import Controller
from blastn import Blastn
from flanking_regions import Flanking_Regions
from predict_pam import Predict_PAM


###############################
#         ARGUMENTS           #
###############################

# required
ap = argparse.ArgumentParser(description='PAMpredict version {}'.format('1.0.0')) #pkg_resources.require("cctyper")[0].version))
ap.add_argument('spacers', help='Input fasta file of CRISPR spacers.')
ap.add_argument('phage_db', help='Directory containing Blastn databases of phagic genomes.')
ap.add_argument('outdir', help='Path to the output directory.')

#optional
ap.add_argument('-t', '--threads', help='Number of parallel processes [%(default)s].', default=1, type=int)
ap.add_argument('--keep_tmp', help='Keep temporary files.', action='store_true')
ap.add_argument('--log_lvl', help='Logging level [%(default)s].', default='INFO', type=str, choices=['DEBUG','INFO','WARNING','ERROR'])
ap.add_argument('--force', help='Overwrites existing results.', action='store_true')
ap.add_argument('-d', '--max_diff', help='Maximum number of differences (gaps + mismatches) allowed between spacers and putative protospacers [%(default)s].', default=4, type=int)
ap.add_argument('-p', '--pam_position', help='PAM position with repsect to spacers, default is [%(default)s] (e.g. for Cas9), can be changed to upstream (e.g. for Cas12)', default='DOWNSTREAM', type=str, choices=['UPSTREAM', 'DOWNSTREAM'])
ap.add_argument('-f', '--format', help='File format of the PAM plot [%(default)s].', type=str, default='png', choices=['png','ps','eps','svg','pdf'])
ap.add_argument('-l', '--pam_length', help='Number of PAM positions use to generate predictions and plot [%(default)s].', default=10, type=int)
ap.add_argument('--no_plot', help='Suppress plot generation.', action='store_true')


###############################
#         RUN PIPLINE         #
###############################

# Initialize
master = Controller(ap.parse_args())

# Blastn
protospacers = Blastn(master)
protospacers.run_blastn()

# Flanking regions realignment
flanking_regions = Flanking_Regions(protospacers)
flanking_regions.run_realignemnt()

# Predict PAM
pam = Predict_PAM(flanking_regions)
pam.predict()
