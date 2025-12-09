#!/bin/bash

# dataset='data/Physiome/aligned_small_odes.json'

# MODELS=(
#     "afp"
#     "feafp"
#     "eplex"
#     "ehc"
#     "proged"
#     "ffx"
#     "pysr"
#     "sindy"
# )

# for i in "${!MODELS[@]}"; do
#     model="${MODELS[$i]}"
#     gpu_id=$((i % 8)) # Assign GPU id from 0 to 7
#     echo "Submitting job for model: $model on GPU $gpu_id"
#     sbatch --time=2-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=gpu odeformer_scripts/run_baseline.sh $dataset $model $gpu_id
# done


sbatch --time=3-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G baselines.sh
