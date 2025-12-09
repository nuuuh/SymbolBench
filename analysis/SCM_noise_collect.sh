#!/bin/bash

# python summary/SCM_baseline_summary.py

models=("gpt-4o-mini" "gpt-oss-20b" "gpt-oss-120b") # "gpt-4o-mini" "gpt-oss-20b" "gpt-oss-120b"
noise=(0.01 0.05 0.1 0.2 0.4)
for model in "${models[@]}"; do
    for n in "${noise[@]}"; do
        folder="analysis_output/"
        eval_type="noise_${n}"
        python summary/SCM_summary.py --model $model --folder $folder --eval_type $eval_type
    done
done









