#!/bin/bash
# Batch-run pipeline jobs: Usage: run_pipeline.sh <start_idx> <end_idx> <batch_size>
export OPENAI_API_KEY="your_openai_api_key_here"

dataset="data/strogatz_extended.json"
model="ChatTS-14B" # ChatTS-14b gpt4o_mini
judge_model="gpt-4.1-nano"
additional_prompt='none'

# # reasoning
exp_name="case_check_base"
reasoning=false
use_context=true
judge_enabled=false
prompt_path="prompts/DE/DE_textual_base.json"

for idx in 5; # 8 27 33 43 49 50 58
do
    python pipeline.py \
        device="auto" \
        model=${model} \
        experiment.symbolic_expression.name=${exp_name} \
        experiment.symbolic_expression.data_path=${dataset} \
        experiment.textual_input=${prompt_path} \
        experiment.symbolic_expression.reasoning=${reasoning} \
        experiment.symbolic_expression.use_context=${use_context} \
        experiment.symbolic_expression.idx=$idx \
        experiment.symbolic_expression.judge.enabled=${judge_enabled} \
        experiment.symbolic_expression.judge.gpt_model=${judge_model} \
        experiment.symbolic_expression.additional_prompt=${additional_prompt}
done


