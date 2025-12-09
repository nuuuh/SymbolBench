#!/bin/bash
# Batch-run pipeline jobs: Usage: run_pipeline.sh <start_idx> <end_idx> <batch_size>
start=${1:-0}
end=${2:-189}
batch_size=${3:-1}

export OPENAI_API_KEY="your_openai_api_key_here"
export TOGETHER_API_KEY="your_together_api_key_here"

 # 'none' 'image'
dataset='data/science_sgc_total.npy' # data/science_sgc_total.npy or data/synthetic_sgc.npy
model="Qwen2.5-14b" # 
judge_model="gpt-4.1-nano"
additional_prompt='none'


# # base
exp_name="SCM_base"
reasoning=false
use_context=false
judge_enabled=false
prompt_path="prompts/SCM/SCM_textual_base.json"



for ((idx=start; idx<=end; idx++)); do
   # Wait until fewer than batch_size jobs are running
   while [ "$(jobs -rp | wc -l)" -ge "$batch_size" ]; do
       sleep 60
   done
   # Skip if results already exist
   if [ ! -f "runs/Qwen2.5-14B/textual/${exp_name}/Expr_${idx}/final_results.npy" ]; then
   echo "Running idx=${idx}..."
       python pipeline.py \
           device="auto" \
           experiment='SCM_exp' \
           dataset='SCM' \
           model=${model} \
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
