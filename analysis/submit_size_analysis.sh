#!/bin/bash

think=true

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu analysis/DE_analysis.sh 1 29 3 'Qwen3-0.6B' $think 'none' 'cuda:0' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu analysis/DE_analysis.sh 1 29 3 'Qwen3-1.7B' $think 'none' 'cuda:1' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu analysis/DE_analysis.sh 1 29 2 'Qwen3-4B' $think 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/DE_analysis.sh 1 29 2 'Qwen3-8B' $think 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/DE_analysis.sh 1 29 1 'Qwen3-14B' $think 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/DE_analysis.sh 1 29 1 'Qwen3-32B' $think 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/DE_analysis.sh 1 29 1 'Qwen3-235B-A22B' $think 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=256G --partition=l40s analysis/DE_analysis.sh 1 29 2 'DeepSeek-V3' False 'none' 'cpu' 
sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=256G --partition=l40s analysis/DE_analysis.sh 1 29 2 'gpt-oss-120b' True 'none' 'cpu' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=256G --partition=l40s analysis/DE_analysis.sh 1 29 2 'gpt-oss-20b' True 'none' 'cpu'
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=256G --partition=l40s analysis/DE_analysis.sh 1 29 2 'DeepSeek-R1-Distill-Llama-70B-free' True 'none' 'cpu' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s analysis/DE_analysis.sh 1 29 2 'Llama-3.3-70B-Instruct-Turbo-Free' False 'none' 'cpu'  
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s analysis/DE_analysis.sh 1 29 2 'exaone-3-5-32b-instruct' False 'none' 'cpu'  

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/DE_analysis.sh 1 29 1 'o4-mini' $think 'none' 'cpu' 

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/DE_analysis.sh 1 29 1 'Qwen3-32B' true 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/DE_analysis.sh 1 29 1 'Qwen3-32B' false 'none' 'auto' 


# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu analysis/BN_analysis.sh 2 'Qwen3-0.6B' true 'none' 'cuda:3' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu analysis/BN_analysis.sh 2 'Qwen3-0.6B' false 'none' 'cuda:2' 

#-----ignored
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu analysis/BN_analysis.sh 3 'Qwen3-1.7B' $think 'none' 'cuda:1' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu analysis/BN_analysis.sh 2 'Qwen3-4B' $think 'none' 'cuda:2' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/BN_analysis.sh 2 'Qwen3-8B' $think 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/BN_analysis.sh 2 'Qwen3-14B' $think 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/BN_analysis.sh 1 'Qwen3-32B' $think 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/BN_analysis.sh 1 'Qwen3-235B-A22B' $think 'none' 'auto' 
#-----

# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu analysis/SCM_analysis.sh 1 'Qwen3-1.7B' true 'none' 'cuda:3' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu analysis/SCM_analysis.sh 1 'Qwen3-1.7B' false 'none' 'cuda:4' 

#-----ignored
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu analysis/SCM_analysis.sh 3 'Qwen3-1.7B' $think 'none' 'cuda:1' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=64G --partition=gpu analysis/SCM_analysis.sh 2 'Qwen3-4B' $think 'none' 'cuda:2' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/SCM_analysis.sh 2 'Qwen3-8B' $think 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/SCM_analysis.sh 2 'Qwen3-14B' $think 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/SCM_analysis.sh 1 'Qwen3-32B' $think 'none' 'auto' 
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=gpu analysis/SCM_analysis.sh 1 'Qwen3-235B-A22B' $think 'none' 'auto' 
#-----