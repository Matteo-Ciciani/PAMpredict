import os
import sys
import logging
from Bio import SeqIO

class Controller:
    def __init__(self, args):
        self.spacers_fasta = args.spacers
        self.phage_db = args.phage_db
        self.outdir = args.outdir
        self.threads = args.threads
        self.log_lvl = args.log_lvl
        self.output = args.outdir
        self.max_diff = args.max_diff
        self.format = args.format
        self.no_plot = args.no_plot
        self.keep_tmp = args.keep_tmp
        self.pam_len = args.pam_length
        self.pam_position = args.pam_position.lower()
        self.force = args.force
        
        # Logger
        logging.basicConfig(format='\033[36m'+'[%(asctime)s] %(levelname)s:'+'\033[0m'+' %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=self.log_lvl)
        logging.info('Running PAMpredict version {}'.format('1.0.1'))
        
        # Check arguments
        self.check_input()
        self.check_output()
        self.check_db()
        self.check_optional()
    
        # Write arguments
        da = vars(args)
        f = open(os.path.join(self.outdir, 'arguments.tsv'), 'w')
        for k, v in da.items():
            f.write('{}:\t{}\n'.format(k, v))
        f.close()
    
    def check_input(self):
        if os.path.isfile(self.spacers_fasta):
            self.check_fasta()
        else:
            logging.error('Could not find input file')
            sys.exit(1)
    
    def check_fasta(self):
        
        # Get spacers
        with open(self.spacers_fasta, 'r') as handle:
            self.spacer_dict = {}
            for fa in SeqIO.parse(handle, 'fasta'):
                if fa.id in self.spacer_dict:
                    logging.error('Duplicate fasta headers detected')
                    sys.exit(1)
                if not all([x in ['A','T','C','G', 'N'] for x in fa.seq.upper()]):
                    logging.error('Spacer sequence ' + fa.id + ' contains invalid bases!')
                    sys.exit(1)
                self.spacer_dict[fa.id] = str(fa.seq.upper())
            
        # Check for numeric headers
        self.num_headers = False
        for i in self.spacer_dict.keys():
            try:
                dump = float(i)
                self.num_headers = True
            except:
                pass
        
        if self.num_headers:
            logging.warning('Numeric fasta headers detected. A prefix is added to the names')
            new_fasta = open(os.path.join(self.out, 'fixed_input.fna'), 'w')
            subprocess.run(['sed', 's/^>/>Spacer/', self.fasta], stdout = new_fasta)
            new_fasta.close()
            self.fasta = os.path.join(self.out, 'fixed_input.fna')
            self.seq_dict = {'Spacer'+str(key): val for key, val in self.seq_dict.items()}
    
    def check_output(self):
        try:
            os.mkdir(self.outdir)
        except FileExistsError:
            if not self.force:
                logging.error('Directory '+self.outdir+' already exists! Use --force to overwrite.')
                sys.exit(1)
            
    def check_optional(self):
        if self.threads < 1:
            logging.error('Number of threads must be at least 1!')
            sys.exit(1)
        if self.max_diff < 0:
            logging.error('Number of differences between spacers and putative protospacers cannot be less than 0!')
            sys.exit(1)
        if self.pam_len < 3 or self.pam_len > 20:
            logging.error('Invalid PAM length!')
            sys.exit(1)
            
    def check_db(self):
        if os.path.isdir(self.phage_db):
            self.all_db_files = os.listdir(self.phage_db)
            self.db_fasta = [x for x in self.all_db_files if x.endswith(
                '.fna') or x.endswith('.fa') or x.endswith('.fasta')]
            if len(self.db_fasta) < 1:
                logging.error('No fasta file found in ' + self.phage_db)
                sys.exit(1)
            #for db in self.db_fasta:
            #    for ext in ['.nhr','.nin','.nog','.nsd','.nsi','.nsq']:
            #        if db + ext not in self.all_db_files:
            #            logging.error('Blast database file not found: ' + db + ext)
            #            sys.exit(1)
        else:
            logging.error('Directory not found: ' + self.phage_db)
            sys.exit(1)
