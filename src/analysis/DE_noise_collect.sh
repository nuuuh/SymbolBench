#!/bin/bash
models=("gpt-oss-20b" "gpt-oss-120b") #"gpt-4o-mini" "gpt-oss-20b" "gpt-oss-120b"
modality='textual'
noise_level=(0.01 0.05 0.1 0.2 0.4)
dataset='strogatz'


for model in "${models[@]}"; do
    for noise in "${noise_level[@]}"; do
        folder="analysis_output/${model}/textual/DE_noise_${noise}/"
        et="noise_${noise}"
        python summary/DE_LLM_preprocess.py --model $model --dataset $dataset --folder $folder
        
        settings=(0 1)
        for ood in "${settings[@]}"; do
            python evaluate_DEs.py --model $model --dataset $dataset --eval_type $et --OOD $ood --LLM 1 --modality $modality --folder $folder   
        done
    done
done









