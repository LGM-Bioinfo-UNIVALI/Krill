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

default_threads = len(os.sched_getaffinity(0))

cprint.ok('#################\n### Krill 1.0 ###\n#################')
cprint.ok('\nPlease cite us!\n')

parser = argparse.ArgumentParser()
parser.add_argument('-noprep','--do_not_prepare_fasta_files',help='Rename fasta files, its headers and store changes in a CSV file for control [DEFAULT: TRUE]',action='store_true')
parser.add_argument('-x','--fasta_extension',help='Fasta files extension [DEFAULT: fasta]',type=str,default='fasta')
parser.add_argument('-t','--threads',help='Trheads to use in analysis [DEFAULT: {}]'.format(default_threads),type=int,default=default_threads)
parser.add_argument('--citation',help='Shows how to cite us',action='store_true')
parser.add_argument('PATH',help='Working path with fasta files',type=str)
args = parser.parse_args()

dbs = [os.path.join(args.PATH,x) for x in os.listdir(args.PATH) if os.path.isdir(os.path.join(args.PATH,x))]
fastas=list(pathlib.Path(args.PATH).rglob('*.{}'.format(args.fasta_extension)))
total_size=int(sum([os.stat(x).st_size for x in fastas])/(10**6))



with tqdm(desc='Total Databases Analysis progress in MB',total=total_size,unit='MB',colour='blue',position=0,leave=False) as dbar:
  for db in dbs:

      os.chdir(db)

      if not args.do_not_prepare_fasta_files:
          prepareFastas.run(db,args.fasta_extension)

      os.makedirs(os.path.join(db,'AntiSMASH'),exist_ok=True)
      os.makedirs(os.path.join(db,'ARTS'),exist_ok=True)

      db_fastas = list(pathlib.Path(db).glob('*.{}'.format(args.fasta_extension)))

      
      with tqdm(desc='{db} BGCs Analysis'.format(db=os.path.basename(db)),total=len(db_fastas),unit='file',colour='blue',position=1,leave=False) as pbar:
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
          result_futures = [executor.submit(run_AntiSMASH_and_ARTS.run, fasta_file, args.fasta_extension, db) for fasta_file in db_fastas]
          for future in as_completed(result_futures):
            try:
              pbar.update(1)
              dbar.update(int(future.result()))
            except Exception as e:
              cprint.fatal(e,interrupt=False)
      

      result_futures = [run_AntiSMASH_and_ARTS.run(fasta_file, args.fasta_extension, db) for fasta_file in db_fastas]
      for future in result_futures:
        try:
          pbar.update(1)
          dbar.update(int(future))
        except Exception as e:
          cprint.fatal(e,interrupt=False)
          
      cprint.info('# Extracting AntiSMASH and ARTS results...')
      AntiSMASH_extractor.prepare_extraction(db)
      ARTS_extractor.prepare_extraction(db)
      AntiSMASH_extractor.protocluster_parse(db,args.threads)
      AntiSMASH_extractor.get_clusters_blast(db,args.threads)
      ARTS_extractor.ARTS_overview(db)
      ARTS_extractor.ARTS_Results_Extraction(db)
      #defs = [ARTS_extractor.ARTS_overview,
      #       ARTS_extractor.ARTS_Results_Extraction]
      
      #with ThreadPoolExecutor(max_workers=args.threads) as executor:
      #     for func in defs:
      #          executor.submit(func, db)
      

      cprint.info('# Screening AntiSMASH and ARTS results...')
      get_screening_results.run(db)

      cprint.info('# Getting fasta info (Bases and ORFs count)...')
      count_bases_and_ORFs.run(db,args.fasta_extension,args.threads)

      os.chdir(args.PATH)

cprint.info('# Getting metainfo (Databases MBases, ORFs, BGCs/Mbases and BGCs/ORFs)...')
get_dbs_metainfo.get(args.PATH,args.fasta_extension)

cprint.ok('\nAll done. Any questions please contact: saulobdasilva@gmail.com\nCheers!')
