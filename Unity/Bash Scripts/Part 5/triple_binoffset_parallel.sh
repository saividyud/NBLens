#!/usr/bin/env bash

#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=triple_lens_array_5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-124
#SBATCH --gres=gpu:1

#SBATCH --output="./Unity/Output Logs/Part_5/triple_%a_output.txt"

# Commands to run
module load mamba
module load cuda/12.9.0-vvrwkrx
export LD_LIBRARY_PATH=$CUDA_HOME/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
mamba activate .venv

# Defining parameters
ss1=(0.2 0.6 1.0 1/0.6 1/0.2)
ss2=(0.2 0.6 1.0 1/0.6 1/0.2)
alphas=(0 45 90 135 180)

# Compute indices (5×5×5 grid: s1 varies slowest, alpha varies fastest)
s1_index=$((SLURM_ARRAY_TASK_ID / 25))
s2_index=$(((SLURM_ARRAY_TASK_ID % 25) / 5))
alpha_index=$((SLURM_ARRAY_TASK_ID % 5))

# Extract parameters
s1=${ss1[$s1_index]}
s2=${ss2[$s2_index]}
alpha=${alphas[$alpha_index]}

echo "Running $s1 separation and $s2 separation and $alpha angle"

python "./Unity/Python Scripts/Part 5/triple_binoffset_parallel.py" -s1 $s1 -s2 $s2 -alpha $alpha