import subprocess, os, datetime

def run(fasta_file, db):
    project_dir = os.path.dirname(os.path.abspath(__file__))
    name=str(os.path.basename(fasta_file)).rsplit('.',1)[0]
    time=datetime.datetime.now().strftime('%Y-%m-%d %X')
    
    # AntiSMASH
    if not os.path.isdir(os.path.join(db,'AntiSMASH',name)):
        pass
        subprocess.run(
            'bash {project_dir}/krill_run_antismash.sh {file} {thr} AntiSMASH AntiSMASH/{file_name}'.format(project_dir=project_dir, file=fasta_file, thr=1, file_name=name),
            shell=True,
            executable='/bin/bash'
        )

        # subprocess.run('''bash run_antismash {file} AntiSMASH \
        #                             --cc-mibig --rre --tigrfam \
        #                             --skip-zip-file \
        #                             --cb-general \
        #                             --cb-knownclusters \
        #                             --cb-subclusters \
        #                             --asf \
        #                             --pfam2go \
        #                             --smcog-trees \
        #                             --genefinding-tool prodigal-m \
        #                             --cpus {thr} 
        #                '''.format(
        #                              file=fasta_file,
        #                              thr=1
        #                              ),
        #                           shell=True,
        #                           executable='/bin/bash')
    else:
        pass
    

    #ARTS
    if not os.path.isdir(os.path.join(db,'ARTS',name)):
        hmm_ref=os.path.join(os.path.dirname(__file__),'db','bacterial_metagenomic_resist_models.hmm')
        res_dir=os.path.join(db,'ARTS',name)
        os.makedirs(res_dir,exist_ok=True)
  
        thr=1
        gbk=os.path.join(db,'AntiSMASH',name,name+'.gbk')

        subprocess.run(
            'bash {project_dir}/krill_run_arts.sh {input} {hmm_ref} {res_dir} {thr}'.format(project_dir=project_dir, input=gbk, hmm_ref=hmm_ref, res_dir=res_dir, thr=thr),
            shell=True,
            executable='/bin/bash',
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT
        )

        # subprocess.run('''source ~/miniconda3/etc/profile.d/conda.sh && \
        #                conda activate ARTS && \
        #                python ~/miniconda3/envs/ARTS/arts/artspipeline1.py \
        #                -org METAGENOME \
        #                -khmms {hmm_ref} \
        #                -t E3 \
        #                -rd {res_dir} \
        #                -cpu {thr} \
        #                {input} \
        #                ~/miniconda3/envs/ARTS/arts/reference/metagenome/
        #                '''.format(
        #                              hmm_ref=hmm_ref,
        #                              res_dir=res_dir,
        #                              thr=thr,
        #                           input=gbk
        #                          ),
        #                          shell=True,
        #                          executable='/bin/bash',
        #                          stdout=subprocess.DEVNULL,
        #                          stderr=subprocess.STDOUT)
    #else:
    #  pass
    return int(os.stat(fasta_file).st_size/(10**6))
