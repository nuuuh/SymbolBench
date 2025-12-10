#!/bin/bash

# python summary/SCM_baseline_summary.py

models=("Llama-3.2-3B" ) # "Qwen2.5-14B"  "gpt-4o-mini" "Llama-3.2-3B" "Mathstral-7B-v0.1"
eval_types=( "context" "reasoning") # "naive" "base" "context" "reasoning"

for model in "${models[@]}"; do
    for eval_type in "${eval_types[@]}"; do
        python summary/SCM_summary.py --model $model --eval_type $eval_type
    done
done









