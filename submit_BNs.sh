#!/bin/bash
# model="Qwen2.5-14B"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image
# partition='l4-8-gm192-c192-m768'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_context.sh 1 65 1 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 cripts/BN/BN_reasoning.sh 1 65 1 $model $additional_prompt 'cuda'


# model="ChatTS-14b"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='time_series' # none, time_series, image
# partition='l4-8-gm192-c192-m768'


# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_context.sh 1 65 1 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 cripts/BN/BN_reasoning.sh 1 65 1 $model $additional_prompt 'cuda'


# model="Llama-3.2-3B"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image
# partition='l4-8-gm192-c192-m768'

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_context.sh 1 65 1 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_reasoning.sh 1 65 1 $model $additional_prompt 'cuda'

# model="Mathstral-7B-v0.1"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image
# partition='l4-8-gm192-c192-m768'

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_context.sh 1 65 1 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_reasoning.sh 1 65 1 $model $additional_prompt 'cuda'


# model="gpt-4o-mini"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image
# partition='l40s'

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_context.sh 1 65 1 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_reasoning.sh 1 65 1 $model $additional_prompt 'cuda'



# model="gpt-oss-120b"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image
# partition='l40s'

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_context.sh 1 65 1 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_reasoning.sh 1 65 1 $model $additional_prompt 'cuda'

# model="DeepSeek-V3"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image
# partition='l40s'

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_context.sh 1 65 1 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/BN/BN_reasoning.sh 1 65 1 $model $additional_prompt 'cuda'


# -------Hybrid Experiments -------

model="gpt-4o-mini"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
additional_prompt='time_series' # none, time_series, image
llm_weights=(0.01 0.1 0.5 1 2)

sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0 scripts/BN/hybrid_GP.sh 1 65 1 $model true 1 $additional_prompt 'cpu'

# for llm_weight in "${llm_weights[@]}"; do
#     sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=h200 scripts/BN/hybrid.sh 1 29 1 $model true $llm_weight $additional_prompt 'cpu'
# done



# model="gpt-4o-mini"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image
# # llm_weights=(0.01 0.1 0.5 1 2)

# sbatch --time=1-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=h200 scripts/BN/hybrid_LLM.sh 1 29 1 $model 1 $additional_prompt 'cpu'