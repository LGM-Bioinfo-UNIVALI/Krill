#!/bin/bash

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

run_conda_func $1 $2 $3 $4
