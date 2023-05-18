#!/bin/bash

run_docker_func () {
    bash run_antismash $1 AntiSMASH \
    --cc-mibig --rre --tigrfam \
    --skip-zip-file \
    --cb-general \
    --cb-knownclusters \
    --cb-subclusters \
    --asf \
    --pfam2go \
    --smcog-trees \
    --genefinding-tool prodigal-m \
    --cpus $2 
}

run_conda_func () {
    source ~/miniconda3/etc/profile.d/conda.sh && \
    conda activate ARTS && \
    python ~/miniconda3/envs/ARTS/arts/artspipeline1.py \
    -org METAGENOME \
    -khmms $2 \
    -t E3 \
    -rd $3 \
    -cpu $4 \
    $1 \
    ~/miniconda3/envs/ARTS/arts/reference/metagenome/
    conda deactivate
}


if ! command -v docker &> /dev/null; then
    run_conda_func $1 $2 $3 $4
    exit
else
    if [[ $(docker images -q arts) ]]; then
        run_docker_func $1 $2 $3
    else
        run_conda_func $1 $2 $3 $4
    fi
fi
