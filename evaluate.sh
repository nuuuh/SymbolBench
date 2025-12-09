#!/bin/bash

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

eval_types=("naive" "base" "context")
datasets=('physiome' 'strogatz') # 'strogatz' 'physiome'

for et in "${eval_types[@]}"; do
    for dataset in "${datasets[@]}"; do
        python summary/DE_LLM_preprocess.py --dataset $dataset --eval_type $et
    done
done


eval_types=("DE_naive" "DE_base" "DE_context")
datasets=('physiome') # 'strogatz' 'physiome'
settings=(0 1)

for et in "${eval_types[@]}"; do
    for dataset in "${datasets[@]}"; do
        for ood in "${settings[@]}"; do
            python evaluate_DEs.py --model LLM_SR --dataset $dataset --eval_type $et --OOD $ood
        done
    done
done


# models=("odeformer" "proged" "pysr" "sindy") # 'strogatz' 'physiome'

# for model in "${models[@]}"; do
#     python summary/DE_summary.py --model $model
# done


eval_types=("DE_naive" "DE_base" "DE_context")
model="LLM_SR"

for et in "${eval_types[@]}"; do
    python summary/DE_summary.py --model $model --eval_type $et
done





