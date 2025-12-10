#!/bin/bash

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate llm



dataset=${1:-'data/strogatz_extended.json'}
model=${2:-'pysr'}
gpu_id=${3:-0}
export CUDA_VISIBLE_DEVICES=$gpu_id

python odeformer_scripts/run_baselines.py \
    --baseline_model="${model}" \
    --dataset="${dataset}"

