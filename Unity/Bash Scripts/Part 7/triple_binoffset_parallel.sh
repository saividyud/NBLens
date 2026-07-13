#!/usr/bin/env bash

#SBATCH --time=02:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=triple_lens_parallel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-59
#SBATCH --gres=gpu:1

#SBATCH --output="./Unity/Output Logs/Part_7/triple_parallel_%a_output.txt"

# Commands to run
module load mamba
module load cuda/12.9.0-vvrwkrx
export LD_LIBRARY_PATH=$CUDA_HOME/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
mamba activate .venv

# Defining parameters (s1 fixed at 1.2)
s1=1.2
ss2=(0.2 0.5 0.8 1/0.8 1/0.5 1/0.2)
q2s=(1e-4 3e-4)
alphas=(0 45 90 135 180)

# Compute indices (6×2×5 grid: s2 varies slowest, alpha varies fastest)
s2_index=$((SLURM_ARRAY_TASK_ID / 10))
q2_index=$(((SLURM_ARRAY_TASK_ID % 10) / 5))
alpha_index=$((SLURM_ARRAY_TASK_ID % 5))

# Extract parameters
s2=${ss2[$s2_index]}
q2=${q2s[$q2_index]}
alpha=${alphas[$alpha_index]}

echo "Running $s1 separation and $s2 separation and $q2 mass ratio and $alpha angle"

python "./Unity/Python Scripts/Part 7/triple_binoffset_parallel.py" -s1 $s1 -s2 $s2 -q2 $q2 -alpha $alpha
