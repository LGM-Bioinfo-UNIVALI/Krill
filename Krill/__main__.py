import sys, os, pathlib, argparse, datetime, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from cprint import *

import prepareFastas
import ARTS_extractor
import AntiSMASH_extractor
import get_screening_results
import count_bases_and_ORFs
import run_AntiSMASH_and_ARTS
import get_dbs_metainfo
import build_charts

default_threads = len(os.sched_getaffinity(0))

cprint.ok('#################\n### Krill 1.0 ###\n#################')
cprint.ok('\nPlease cite us!\n')

parser = argparse.ArgumentParser()
parser.add_argument('-noprep','--do_not_prepare_fasta_files',help='Rename fasta files, its headers and store changes in a CSV file for control [DEFAULT: TRUE]',action='store_true')
parser.add_argument('-t','--threads',help='Threads to use in analysis [DEFAULT: {}]'.format(default_threads),type=int,default=default_threads)
parser.add_argument('--citation',help='Shows how to cite us',action='store_true')
parser.add_argument('PATH',help='Working path with fasta files',type=str)
args = parser.parse_args()



# Get databases (directories) paths if the root directory has them, otherwise work with the root itself
dbs = []
has_multi_db = False
for x in os.listdir(args.PATH):
  if os.path.isdir(os.path.join(args.PATH,x)):
    dbs.append(os.path.join(args.PATH,x))
    has_multi_db = True
  else:
    dbs.append(args.PATH)
    break


# Convert all files extensions (fa, fna, etc) to "fasta"
prepareFastas.convert2fasta(args.PATH)


# Calculate total size of data to be analyzed
fastas=list(pathlib.Path(args.PATH).rglob('*.{}'.format("fasta")))
total_size=int(sum([os.stat(x).st_size for x in fastas])/(10**6))


# Display a visual progress bar for the processing of the whole execution
with tqdm(desc='Total Databases Analysis progress in MB',total=total_size,unit='MB',colour='blue',position=0,leave=False) as dbar:
  for db in dbs:  # Process each database/directory

      os.chdir(db)
    
      if not args.do_not_prepare_fasta_files:
          prepareFastas.run(db,"fasta")  # Rename fasta files, its headers and store changes in a CSV file for control

      os.makedirs(os.path.join(db,'AntiSMASH'),exist_ok=True)  # Create ARTS folder
      os.makedirs(os.path.join(db,'ARTS'),exist_ok=True)  # Create ARTS folder
  
      os.makedirs(os.path.join(db,'ReportOutput'),exist_ok=True)  # Create ARTS folder
    
      db_fastas = list(pathlib.Path(db).glob('*.{}'.format("fasta")))  # List of fasta files

      # Display a visual progress bar for the processing of the current database
      with tqdm(desc='{db} BGCs Analysis'.format(db=os.path.basename(db)),total=len(db_fastas),unit='file',colour='blue',position=1,leave=False) as pbar:
        # Use multiple threads for the execution of the fasta files
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
          # Run AntiSMASH and ARTS
          result_futures = [executor.submit(run_AntiSMASH_and_ARTS.run, fasta_file, db) for fasta_file in db_fastas]
          # Update progress bar
          for future in as_completed(result_futures):
            try:
              pbar.update(1)
              dbar.update(int(future.result()))
            except Exception as e:
              cprint.fatal(e,interrupt=False)

      cprint.info('# Extracting AntiSMASH and ARTS results...')
      AntiSMASH_extractor.prepare_extraction(db)  # Create a folder for the AntiSMASH extraction
      ARTS_extractor.prepare_extraction(db) # Create a folder for the ARTS extraction
      
      AntiSMASH_extractor.protocluster_parse(db, args.threads)  # Parse candidate clusters and coding sequences. output: clusters.tsv
      AntiSMASH_extractor.get_clusters_blast(db, args.threads)
      
      ARTS_extractor.ARTS_overview(db)
      ARTS_extractor.ARTS_Results_Extraction(db)
    

      cprint.info('# Screening AntiSMASH and ARTS results...')
      get_screening_results.run(db)

      cprint.info('# Getting fasta info (Bases and ORFs count)...')
      count_bases_and_ORFs.run(db, "fasta", args.threads)
      
      os.chdir(args.PATH)


os.makedirs(os.path.join(args.PATH,'DBsReportOutput'),exist_ok=True)  # Create ARTS folder

cprint.info('# Getting metainfo (Databases MBases, ORFs, BGCs/Mbases and BGCs/ORFs)...')
get_dbs_metainfo.get(args.PATH, "fasta", has_multi_db)

cprint.info('# Building charts...')
build_charts.build_charts(args.PATH)

cprint.ok('\nAll done. Any questions please contact: saulobdasilva@gmail.com or ellen.junker@edu.univali.br\nCheers!')
