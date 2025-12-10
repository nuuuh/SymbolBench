#!/bin/bash

start=1
end=20
batch_size=3
additional_prompt='none' # 'none' 'time_series' 'image'

models=('gpt-oss-120b' 'gpt-oss-20b') # 'gpt-4o-mini' 'gpt-oss-120b' 'gpt-oss-20b'
noise_level=(0.1 0.2 0.4) # 0.01 0.05 0.1 0.2 0.4
for model in "${models[@]}"; do
    echo "Submitting jobs for model: $model"
    for noise in "${noise_level[@]}"; do
        echo "Submitting jobs for noise level: $noise"
        sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s analysis/noise_DE.sh $start $end $batch_size $model $additional_prompt $noise 'cpu'
        # sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s analysis/noise_SCM.sh $start $end $batch_size $model $additional_prompt $noise 'cpu'
        # sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s analysis/noise_BN.sh $start $end $batch_size $model $additional_prompt $noise 'cpu'
    done
done



# model='gpt-oss-20b'
# additional_prompt='none'
# noise=0.4
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s analysis/noise_SCM.sh 1 20 1 $model $additional_prompt $noise 'cpu'

# model='gpt-oss-120b'
# additional_prompt='none'
# noise=0.2
# sbatch --time=5-24:00:00 --ntasks=1 --cpus-per-task=8 --mem=32G --partition=l40s analysis/noise_SCM.sh 1 20 1 $model $additional_prompt $noise 'cpu'



