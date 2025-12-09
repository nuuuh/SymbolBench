#!/bin/bash

#----- Boolean Network ------

# python evaluate_bn.py --OOD 1
# python evaluate_bn.py --OOD 0

# #----- Structured Causal Graph  ------

# baselines=( 'pcmci' 'lpcmci'  'j-pcmci+'  ) # 'rpcmci'
# datasets=("science_sgc_total") #  "synthetic_sgc"
# for baseline in "${!baselines[@]}"; do
#     for dataset in "${!datasets[@]}"; do
#         echo "Running baseline: ${baselines[$baseline]} on dataset: ${datasets[$dataset]}"
#         python evaluate_scm.py --baseline "${baselines[$baseline]}" --dataset "${datasets[$dataset]}"
#     done
# done

# #----- Scientific Equations ------

# dataset='data/strogatz_extended.json'
# dataset='data/Physiome/solved_small_odes.json'

# MODELS=(
#     # "ffx"
#     # "afp"
#     "pysr"
#     "sindy"
#     "feafp"
#     "eplex"
#     "ehc"
#     "proged"
    
# )

# for i in "${!MODELS[@]}"; do
#     model="${MODELS[$i]}"
#     bash odeformer_scripts/run_baseline.sh $dataset $model 0
# done
