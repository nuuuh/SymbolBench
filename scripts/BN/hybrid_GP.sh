#!/bin/bash
# Batch-run pipeline jobs: Usage: run_pipeline.sh <start_idx> <end_idx> <batch_size>
start=${1:-1}
end=${2:-65}
batch_size=${3:-1}
model=${4:-'Qwen2.5-14B'} # ChatTS-14b
use_llm=${5:-true} # true or false
llm_weight=${6:-1} # Weight for the LLM score in the hybrid model
additional_prompt=${7:-'time_series'} # 'none' 'time_series' 'image'
device=${8:-'auto'} # 'auto' or specific device like 'cuda:0'


export OPENAI_API_KEY="your_openai_api_key_here"
export TOGETHER_API_KEY="your_together_api_key_here"
 # 'none' 'time_series' 'image'
# dataset="data/Physiome/small_odes.json"
dataset="data/BN_small.json"

# context
if $use_llm; then
    exp_name="hybrid_GP_BN_w_llm_${llm_weight}"
else
    exp_name="hybrid_GP_BN_wo_llm"
fi

prompt_path="prompts/BN/BN_hybrid_score.json"


for ((idx=start; idx<=end; idx++)); do
   # Wait until fewer than batch_size jobs are running
   while [ "$(jobs -rp | wc -l)" -ge "$batch_size" ]; do
       sleep 60
   done
   if [ ! -f "hybrid_output/${model}/textual/${exp_name}/Expr_${idx}/final_results.npy" ]; then
        echo "Running idx=${idx}..."
    python hybrid_pipeline.py \
        output_dir="hybrid_output" \
        device=$device \
        model=${model} \
        experiment="BN_exp" \
        dataset='BN' \
        use_llm=$use_llm \
        llm_weight=$llm_weight \
        experiment.symbolic_expression.use_context=true \
        experiment.textual_input=${prompt_path} \
        experiment.symbolic_expression.name=${exp_name} \
        experiment.symbolic_expression.data_path=${dataset} \
        experiment.symbolic_expression.idx=$idx \
        experiment.symbolic_expression.additional_prompt=${additional_prompt} &
    else
        echo "Skipping idx=${idx}, results already exist."
    fi
 done

