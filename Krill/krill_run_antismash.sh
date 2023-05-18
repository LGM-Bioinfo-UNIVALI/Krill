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
    antismash --cc-mibig --rre --tigrfam \
    --skip-zip-file \
    --cb-general \
    --cb-knownclusters \
    --cb-subclusters \
    --asf \
    --pfam2go \
    --smcog-trees \
    --genefinding-tool prodigal-m \
    --cpus $2 \
    --output-dir $3 \
    $1
}


if ! command -v docker &> /dev/null; then
    run_conda_func $1 $2 $4
    exit
else
    if [[ $(docker images -q antismash/standalone) ]]; then
        run_docker_func $1 $2 $3
    else
        run_conda_func $1 $2 $4
    fi
fi
