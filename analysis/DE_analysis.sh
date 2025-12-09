#!/bin/bash
# Batch-run pipeline jobs: Usage: run_pipeline.sh <start_idx> <end_idx> <batch_size>
start=${1:-1}
end=${2:-29}
batch_size=${3:-1}
model=${4:-'Qwen3-0.6B'} # Llama-2-7b  Llama-2-13b-chat-hf
think=${5:-true} # true or false
additional_prompt=${6:-'none'} # 'none' 'time_series' 'image'
device=${7:-'cuda'} # 'auto' or specific device like 'cuda:0'

export OPENAI_API_KEY="your_openai_api_key_here"
export TOGETHER_API_KEY="your_together_api_key_here"
skipped=(2)

 # 'none' 'time_series' 'image'
dataset="data/DE_4_dims.json"
judge_model="gpt-4.1-nano"

# # reasoning
if [[ $think == true ]]; then
    exp_name="DE_size_compute_w_thinking"
else
    exp_name="DE_size_compute_wo_thinking"
fi

reasoning=false
use_context=true
judge_enabled=false
prompt_path="prompts/DE/DE_textual_context.json"

# meta-llama/Llama-2-7b

for ((idx=start; idx<=end; idx++)); do
    if [[ " ${skipped[@]} " =~ " ${idx} " ]]; then
        echo "Skipping idx=${idx}..."
        continue
    fi
   # Wait until fewer than batch_size jobs are running
    while [ "$(jobs -rp | wc -l)" -ge "$batch_size" ]; do
       sleep 60
    done
    if [ ! -f "analysis_output/${model}/textual/${exp_name}/Expr_${idx}/final_results.npy" ]; then
        echo "Running idx=${idx}..."
        python pipeline.py \
            output_dir="analysis_output" \
            device=$device \
            model=${model} \
            enable_thinking=${think} \
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
# Wait for all tasks to complete
wait

