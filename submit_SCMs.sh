#!/bin/bash
# model="Qwen2.5-14B"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image
# partition='a100-8-gm320-c96-m1152'
# # sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_context.sh 0 189 2 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_reasoning.sh 0 40 1 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_reasoning.sh 41 82 1 $model $additional_prompt 'cuda'
# partition='h200-8-gm1128-c192-m2048'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_reasoning.sh 83 123 1 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_reasoning.sh 124 189 1 $model $additional_prompt 'cuda'

# model="ChatTS-14B"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='time_series' # none, time_series, image
# partition='h200'

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_context.sh 0 189 1 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_reasoning.sh 0 189 1 $model $additional_prompt 'cuda'

# model="Llama-3.2-3B"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image
# partition='h200'

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_context.sh 0 189 6 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_reasoning.sh 0 189 6  $model $additional_prompt 'cuda'

# model="Mathstral-7B-v0.1"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image
# partition='h200'

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_context.sh 0 189 3 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_reasoning.sh 0 189 3 $model $additional_prompt 'cuda'


# model="gpt-4o-mini"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
# additional_prompt='none' # none, time_series, image

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G scripts/SCM/SCM_context.sh 0 189 3 $model $additional_prompt 'cuda'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G scripts/SCM/SCM_reasoning.sh 0 189 3 $model $additional_prompt 'cuda'



model="gpt-oss-120b"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
additional_prompt='none' # none, time_series, image
partition='l40s'

sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_context.sh 0 189 3 $model $additional_prompt 'cpu'
sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_reasoning.sh 0 189 3 $model $additional_prompt 'cpu'

model="DeepSeek-V3"  # gpt-4o-mini Qwen2.5-14B ChatTS-14B Llama-3.2-3B Mathstral-7B-v0.1
additional_prompt='none' # none, time_series, image
partition='l40s'

sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_context.sh 0 189 3 $model $additional_prompt 'cpu'
sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=4 --mem=32G --partition=$partition --gpus=1 scripts/SCM/SCM_reasoning.sh 0 189 3 $model $additional_prompt 'cpu'