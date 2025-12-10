#!/bin/bash

# python summary/SCM_baseline_summary.py

models=("gpt-4o-mini" ) # "Qwen2.5-14B"  "gpt-4o-mini" "Llama-3.2-3B" "Mathstral-7B-v0.1"
eval_types=( "noise_0.01" "noise_0.05" "noise_0.1" "noise_0.2" "noise_0.4") # "naive" "base" "context" "reasoning"

for model in "${models[@]}"; do
    for eval_type in "${eval_types[@]}"; do
        python summary/SCM_summary.py --folder analysis_output --model $model --eval_type $eval_type
    done
done









