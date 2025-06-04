#!/bin/bash
#SBATCH --job-name=hmsc-trial
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --gpus-per-node=1
#SBATCH --time=15:00:00
#SBATCH --array=0-3%4

# In your sbatch script, launch in the background
nvidia-smi --query-gpu=timestamp,index,name,utilization.gpu,memory.used,memory.total --format=csv -l 10 > gpu_usage_${SLURM_ARRAY_TASK_ID}.log &
NVIDIA_MONITOR_PID=$!

module restore myhmscstack
# pip install hmsc ...

source params.txt 
SAM=$nSamples
THIN=$thin

input_path="init_file.rds"
output_path="post_file.rds"
output_path=$(printf "post_chain%.2d_file.rds" $SLURM_ARRAY_TASK_ID)

srun python3 -m hmsc.run_gibbs_sampler \
  --input $input_path \
  --output $output_path \
  --samples $SAM \
  --transient $(($SAM*$THIN)) \
  --thin $THIN \
  --verbose 100 \
  --chain $SLURM_ARRAY_TASK_ID

# Kill monitor after job finishes
kill $NVIDIA_MONITOR_PID