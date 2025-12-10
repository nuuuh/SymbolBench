#!/bin/bash
# Batch-run pipeline jobs: Usage: run_pipeline.sh <start_idx> <end_idx> <batch_size>
start=${1:-0}
end=${2:-20}
batch_size=${3:-1}
model=${4:-'gpt-4o-mini'}
additional_prompt=${5:-'none'} # 'none' 'time_series' 'image'
noise_level=${6:-0.01} # noise level to inject into the data
device=${7:-'cuda:0'} # 'auto' or specific device like 'cuda:0'

export OPENAI_API_KEY="your_openai_api_key_here"
export TOGETHER_API_KEY="your_together_api_key_here"
export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7


indices=(42 45 46 47 50 58 60 61 55 1 2 5 6 11 14 21 23 25 27 29) # 20 cases

 # 'none' 'time_series' 'image'
dataset="data/BN_small.json"
judge_model="gpt-4.1-nano"

exp_name="BN_noise_${noise_level}"

reasoning=false
use_context=true
judge_enabled=false
prompt_path="prompts/BN/BN_textual_context.json"

for idx in "${indices[@]}"; do
   # Wait until fewer than batch_size jobs are running
    while [ "$(jobs -rp | wc -l)" -ge "$batch_size" ]; do
       sleep 60
    done
    if [ ! -f "analysis_output/${model}/textual/${exp_name}/Expr_${idx}/final_results.npy" ]; then
        echo "Running idx=${idx}..."
        python pipeline.py \
            experiment='BN_exp' \
            dataset='BN' \
            output_dir="analysis_output" \
            device=$device \
            model=${model} \
            noise_level=$noise_level \
            experiment.symbolic_expression.name=${exp_name} \
            experiment.symbolic_expression.data_path=${dataset} \
            experiment.textual_input=${prompt_path} \
            experiment.symbolic_expression.reasoning=${reasoning} \
            experiment.symbolic_expression.use_context=${use_context} \
            experiment.symbolic_expression.idx=$idx \
            experiment.symbolic_expression.judge.enabled=${judge_enabled} \
            experiment.symbolic_expression.judge.gpt_model=${judge_model} \
            experiment.symbolic_expression.additional_prompt=${additional_prompt} \
            experiment.symbolic_expression.iterations=50 &
    else
        echo "Skipping idx=${idx}, result already exists."
    fi
 done
# Wait for all tasks to complete
wait

