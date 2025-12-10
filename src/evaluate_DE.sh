#!/bin/bash

# baselines
# models=( "odeformer" "proged" "pysr" "sindy" ) #  
# datasets=('strogatz_extended' 'physiome') # 
# settings=(0 1)

# for model in "${models[@]}"; do
#     for dataset in "${datasets[@]}"; do
#         for ood in "${settings[@]}"; do
#             python evaluate_DEs.py --model $model --dataset $dataset --eval_type "" --OOD $ood
#         done
#     done
# done

# models=("odeformer" "proged" "pysr" "sindy") # 'strogatz' 'physiome'

# for model in "${models[@]}"; do
#     python summary/DE_summary.py --model $model
# done


# LLMs

models=("Llama-3.2-3B" ) # "Qwen2.5-14B" "ChatTS-14B" "gpt-4o-mini" "Mathstral-7B-v0.1" "DeepSeek-V3"
modality='textual'


for model in "${models[@]}"; do
    eval_types=("context" "reasoning") # "naive" "base" 
    datasets=( 'strogatz' 'physiome') # 'strogatz' 'physiome'
    
    for et in "${eval_types[@]}"; do
        for dataset in "${datasets[@]}"; do

            python summary/DE_LLM_preprocess.py --model $model --dataset $dataset --eval_type $et --modality $modality
            
            settings=(0 1)
            for ood in "${settings[@]}"; do
                python evaluate_DEs.py --model $model --dataset $dataset --eval_type $et --OOD $ood --LLM 1 --modality $modality
            done
        done

        python summary/DE_summary.py --model $model --eval_type $et --modality $modality

    done
done









