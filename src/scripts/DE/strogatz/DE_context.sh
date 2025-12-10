#!/bin/bash
# Batch-run pipeline jobs: Usage: run_pipeline.sh <start_idx> <end_idx> <batch_size>
start=${1:-1}
end=${2:-63}
batch_size=${3:-1}
model=${4:-'Qwen2.5-14b'} # ChatTS-14b
additional_prompt=${5:-'none'} # 'none' 'time_series' 'image'
device=${6:-'auto'} # 'auto' or specific device like 'cuda:0'


export OPENAI_API_KEY="your_openai_api_key_here"
export TOGETHER_API_KEY="your_together_api_key_here"

 # 'none' 'time_series' 'image'
# dataset="data/Physiome/small_odes.json"
dataset="data/strogatz_extended.json"
judge_model="gpt-4.1-nano"


# context
exp_name="context_strogatz"
reasoning=false
use_context=true
judge_enabled=false
prompt_path="prompts/DE/DE_textual_context.json"


for ((idx=start; idx<=end; idx++)); do
   # Wait until fewer than batch_size jobs are running
    while [ "$(jobs -rp | wc -l)" -ge "$batch_size" ]; do
       sleep 60
    done
    if [ ! -f "runs/${model}/textual/${exp_name}/Expr_${idx}/final_results.npy" ]; then
        echo "Running idx=${idx}..."
        python pipeline.py \
            device=$device \
            model=${model} \
            experiment.symbolic_expression.name=${exp_name} \
            experiment.symbolic_expression.data_path=${dataset} \
            experiment.textual_input=${prompt_path} \
            experiment.symbolic_expression.reasoning=${reasoning} \
            experiment.symbolic_expression.use_context=${use_context} \
            experiment.symbolic_expression.idx=$idx \
            experiment.symbolic_expression.judge.enabled=${judge_enabled} \
            experiment.symbolic_expression.judge.gpt_model=${judge_model} \
            experiment.symbolic_expression.additional_prompt=${additional_prompt} &
    else
       echo "Skipping idx=${idx}, result already exists."
    fi
 done


