#!/bin/bash
# Batch-run pipeline jobs: Usage: run_pipeline.sh <start_idx> <end_idx> <batch_size>
idx=${1:-0}
batch_size=${2:-1}
model=${3:-'Qwen2.5-14B'} # ChatTS-14b
llm_weight=${4:-1} # Weight for the LLM score in the hybrid model
additional_prompt=${5:-'time_series'} # 'none' 'time_series' 'image'
device=${6:-'auto'} # 'auto' or specific device like 'cuda:0'


export OPENAI_API_KEY="your_openai_api_key_here"
export TOGETHER_API_KEY="your_together_api_key_here"

 # 'none' 'time_series' 'image'



if [ "$idx" -eq 0 ]; then
    idx_list=(6 18 19 20 21 26 
              28 29 30 32 34 41 52 62 64 65 75 80 
              81 83 84 89 90 92 93 97 98 100 101 105
              106 107 112 113 127 139 140 141 146 147 148
              153 157 176 177 190 200 202 228 255 
              260 263 266 279 280 292 296 297 299 313 318 
              320 321 322 326 328 331 332 333 334 335 336 339 341 
              348 349 352 357 358 359 362 387 400 425 432 441 
              444 445 456 457 486 489)
else
    idx_list=($idx)
    batch_size=1
fi

exp_name="hybrid_DE_physiome_w_llm"
prompt_path="prompts/DE/DE_hybrid.json"
dataset="data/Physiome/small_odes.json"

reasoning=false
use_context=true
judge_enabled=false


for idx in "${idx_list[@]}"; do
   # Wait until fewer than batch_size jobs are running
    while [ "$(jobs -rp | wc -l)" -ge "$batch_size" ]; do
       sleep 60
    done
    if [ ! -f "hybrid_output/${model}/textual/${exp_name}/Expr_${idx}/final_results.npy" ]; then
        echo "Running idx=${idx}..."
        python pipeline.py \
            output_dir="hybrid_output" \
            device=$device \
            model=${model} \
            experiment="DE_exp" \
            dataset='DE' \
            use_llm=$use_llm \
            llm_weight=$llm_weight \
            hybrid=true \
            experiment.symbolic_expression.name=${exp_name} \
            experiment.symbolic_expression.data_path=${dataset} \
            experiment.textual_input=${prompt_path} \
            experiment.visual_input=${prompt_path} \
            experiment.symbolic_expression.reasoning=${reasoning} \
            experiment.symbolic_expression.use_context=${use_context} \
            experiment.symbolic_expression.idx=$idx \
            experiment.symbolic_expression.judge.enabled=${judge_enabled} \
            experiment.symbolic_expression.judge.gpt_model=${judge_model} \
            experiment.symbolic_expression.additional_prompt=${additional_prompt} &
    else
        echo "Skipping idx=${idx}, results already exist."
    fi
done

