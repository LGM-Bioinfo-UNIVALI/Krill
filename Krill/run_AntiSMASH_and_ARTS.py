import subprocess, os, datetime

def run(fasta_file,ext,db):
    name=str(os.path.basename(fasta_file)).rsplit('.',1)[0]
    time=datetime.datetime.now().strftime('%Y-%m-%d %X')

    #AntiSMASH
    # if not os.path.isdir(os.path.join(db,'AntiSMASH',name)):
    #     subprocess.run('''source ~/miniconda3/etc/profile.d/conda.sh && \
    #                    conda activate antismash6 && \
    #                    antismash --output-dir AntiSMASH/{dest} \
    #                                 --cc-mibig --rre --tigrfam \
    #                                 --skip-zip-file \
    #                                 --cb-general \
    #                                 --cb-knownclusters \
    #                                 --cb-subclusters \
    #                                 --asf \
    #                                 --pfam2go \
    #                                 --smcog-trees \
    #                                 --genefinding-tool prodigal-m \
    #                                 --cpus {thr} \
    #                                 {file}
    #                    '''.format(
    #                                  dest=name,
    #                                  file=fasta_file,
    #                                  thr=1
    #                                  ),
    #                               shell=True,
    #                               executable='/bin/bash')
    # #else:
    # #  pass

    # AntiSMASH
    if not os.path.isdir(os.path.join(db,'AntiSMASH',name)):
        subprocess.run('''bash run_antismash {file} AntiSMASH \
                                    --cc-mibig --rre --tigrfam \
                                    --skip-zip-file \
                                    --cb-general \
                                    --cb-knownclusters \
                                    --cb-subclusters \
                                    --asf \
                                    --pfam2go \
                                    --smcog-trees \
                                    --genefinding-tool prodigal-m \
                                    --cpus {thr} 
                       '''.format(
                                     file=fasta_file,
                                     thr=1
                                     ),
                                  shell=True,
                                  executable='/bin/bash')


    else:
     pass
    

    #ARTS
    if not os.path.isdir(os.path.join(db,'ARTS',name)):
     
        hmm_ref=os.path.join(os.path.dirname(__file__),'db','bacterial_metagenomic_resist_models.hmm')
        res_dir=os.path.join(db,'ARTS',name)
        os.makedirs(res_dir,exist_ok=True)
  
        thr=1
        gbk=os.path.join(db,'AntiSMASH',name,name+'.gbk')

        subprocess.run('''source ~/miniconda3/etc/profile.d/conda.sh && \
                       conda activate ARTS && \
                       python ~/miniconda3/envs/ARTS/arts/artspipeline1.py \
                       -org METAGENOME \
                       -khmms {hmm_ref} \
                       -t E3 \
                       -rd {res_dir} \
                       -cpu {thr} \
                       {input} \
                       ~/miniconda3/envs/ARTS/arts/reference/metagenome/
                       '''.format(
                                     hmm_ref=hmm_ref,
                                     res_dir=res_dir,
                                     thr=thr,
                                  input=gbk
                                 ),
                                 shell=True,
                                 executable='/bin/bash',
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.STDOUT)
    #else:
    #  pass
    return int(os.stat(fasta_file).st_size/(10**6))





#source ~/miniconda3/etc/profile.d/conda.sh && conda activate ARTS && python ~/miniconda3/envs/ARTS/arts/artspipeline1.py -org METAGENOME -khmms Krill/db/bacterial_metagenomic_resist_models.hmm -t E3 -rd /ssdRaid/iscsi_data/saulo/ellen/Deception_Island/ARTS/000000001 -cpu 1 /ssdRaid/iscsi_data/saulo/ellen/Deception_Island/AntiSMASH/000000001/000000001.gbk ~/miniconda3/envs/ARTS/arts/reference/metagenome/