# #!/bin/bash

export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7

model="ChatTS-14B"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
additional_prompt='time_series' # none, time_series, image

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_naive.sh 1 63 1 $model 'none' 'cuda:0'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=3-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0   scripts/DE/strogatz/DE_base.sh 1 63 2 $model 'none' 'cuda:0'
sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu --gres=gpu:1 --gpu-bind=map_gpu:5   scripts/DE/strogatz/DE_context.sh 1 63 1 $model $additional_prompt
sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu --gres=gpu:1 --gpu-bind=map_gpu:6   scripts/DE/strogatz/DE_reasoning.sh 1 63 1 $model $additional_prompt


# sbatch --time=6-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:1   scripts/DE/physiome/DE_naive.sh 0 1 $model 'none' 'cuda:0'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:2   scripts/DE/physiome/DE_base.sh 0 1 $model 'none' 'cuda:0'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_context.sh 0 1 $model $additional_prompt 'cuda:0'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_reasoning.sh 0 1 $model $additional_prompt 'cuda:0'


# model="Llama-3.2-3B"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:3   scripts/DE/strogatz/DE_naive.sh 1 63 1 $model 'none' 'cuda:0'
# sbatch --time=3-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:3   scripts/DE/strogatz/DE_base.sh 1 63 1 $model 'none' 'cuda:0'
# # CUDA_VISIBLE_DEVICES=0 sbatch --time=3-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0    scripts/DE/strogatz/DE_context.sh 1 31 1 $model $additional_prompt 'cuda:0'
# # CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0    scripts/DE/strogatz/DE_reasoning.sh 1 31 1 $model $additional_prompt 'cuda:0'


# sbatch --time=6-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:4   scripts/DE/physiome/DE_naive.sh 0 1 $model 'none' 'cuda:0'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:4    scripts/DE/physiome/DE_base.sh 0 1 $model 'none' 'cuda:0'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_context.sh 0 1 $model $additional_prompt 'cuda:0'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_reasoning.sh 0 1 $model $additional_prompt 'cuda:0'



# model="Mathstral-7B-v0.1"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:5  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_naive.sh 1 63 1 $model 'none' 'cuda:0'
# sbatch --time=3-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:5  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_base.sh 1 63 1 $model 'none' 'cuda:0'
# # CUDA_VISIBLE_DEVICES=0 sbatch --time=3-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_context.sh 1 31 1 $model $additional_prompt 'cuda:0'
# # CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_reasoning.sh 1 31 1 $model $additional_prompt 'cuda:0'


# sbatch --time=6-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:6  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_naive.sh 0 1 $model 'none' 'cuda:0'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:6  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_base.sh 0 1 $model 'none' 'cuda:0'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_context.sh 0 1 $model $additional_prompt 'cuda:0'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_reasoning.sh 0 1 $model $additional_prompt 'cuda:0'


# model="gpt-4o-mini"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image

# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G scripts/BN/BN_context.sh 1 65 4 $model $additional_prompt 'cpu'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G scripts/BN/BN_reasoning.sh 1 65 4 $model $additional_prompt 'cpu'

# model="gpt-4o-mini"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='image' # none, time_series, image

# # CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G scripts/DE/strogatz/DE_context.sh 1 63 4 $model $additional_prompt 'cpu'
# # CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G scripts/DE/strogatz/DE_reasoning.sh 1 63 4 $model $additional_prompt 'cpu'

# # CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G scripts/DE/physiome/DE_context.sh 0 2 $model $additional_prompt 'cpu'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G scripts/DE/physiome/DE_reasoning.sh 0 4 $model $additional_prompt 'cpu'

model="Qwen2.5-14B"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
additional_prompt='none' # none, time_series, image

# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=turing --gres=gpu:1 --gpu-bind=map_gpu:0 scripts/DE/strogatz/DE_context.sh 1 63 2 $model $additional_prompt 'cuda:0'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=turing --gres=gpu:1 --gpu-bind=map_gpu:0 scripts/DE/strogatz/DE_reasoning.sh 1 63 2 $model $additional_prompt 'cuda:0'

sbatch --time=6-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_naive.sh 0 1 $model $additional_prompt
sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu --gres=gpu:1 --gpu-bind=map_gpu:1   scripts/DE/physiome/DE_base.sh 0 1 $model $additional_prompt
sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu --gres=gpu:1 --gpu-bind=map_gpu:2 scripts/DE/physiome/DE_context.sh 0 1 $model $additional_prompt 
sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu --gres=gpu:1 --gpu-bind=map_gpu:3 scripts/DE/physiome/DE_reasoning.sh 0 2 $model $additional_prompt



# model="DeepSeek-V3"  
# additional_prompt='none' # none, time_series, image

# CUDA_VISIBLE_DEVICES=0 sbatch --time=3-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_context.sh 1 31 1 $model $additional_prompt 'cpu'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_reasoning.sh 1 31 1 $model $additional_prompt 'cpu'

# CUDA_VISIBLE_DEVICES=0 sbatch --time=3-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_context.sh 32 63 1 $model $additional_prompt 'cpu'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_reasoning.sh 32 63 1 $model $additional_prompt 'cpu'

# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_context.sh 0 1 $model $additional_prompt 'cpu'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_reasoning.sh 0 1 $model $additional_prompt 'cpu'




# model="gpt-oss-120b"  
# additional_prompt='none' # none, time_series, image

# CUDA_VISIBLE_DEVICES=0 sbatch --time=3-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_context.sh 1 31 1 $model $additional_prompt 'cpu'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_reasoning.sh 1 31 1 $model $additional_prompt 'cpu'

# CUDA_VISIBLE_DEVICES=0 sbatch --time=3-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_context.sh 32 63 1 $model $additional_prompt 'cpu'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/strogatz/DE_reasoning.sh 32 63 1 $model $additional_prompt 'cpu'

# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_context.sh 0 3 $model $additional_prompt 'cpu'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0  scripts/DE/physiome/DE_reasoning.sh 0 3 $model $additional_prompt 'cpu'


# -------Hybrid Experiments -------


# model="gpt-4o-mini"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='time_series' # none, time_series, image
# llm_weights=(0.01 0.1 0.5 1 2)

# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0 scripts/DE/hybrid_GP.sh 1 23 1 $model false 1 $additional_prompt 'cpu'
# CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0 scripts/DE/hybrid_GP.sh 1 23 1 $model true 0.1 $additional_prompt 'cpu'
# for llm_weight in "${llm_weights[@]}"; do
#     CUDA_VISIBLE_DEVICES=0 sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0  --gres=gpu:1 --gpu-bind=map_gpu:0 scripts/DE/hybrid.sh 1 28 1 $model true $llm_weight $additional_prompt 'cpu'
# done


# model="gpt-4o-mini"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image
# llm_weights=(0.01 0.1 0.5 1 2)

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=128G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0 scripts/DE/hybrid_LLM.sh 1 63 3 $model 1 $additional_prompt 'cpu'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=128G --partition=l40s --gres=gpu:1 --gpu-bind=map_gpu:0 scripts/DE/hybrid_LLM_2.sh 0 3 $model 1 $additional_prompt 'cpu'



