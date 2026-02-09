#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=binary_animator_4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-118

#SBATCH --output="./Unity/Output Logs/Analysis/animator_%a_output.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
ss=(0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1/0.9 1/0.8 1/0.7 1/0.6 1/0.5 1/0.4 1/0.3 1/0.2)
qs=(1e-6 3e-6 1e-5 3e-5 1e-4 3e-4 1e-3)

# Compute indices
s_index=$((SLURM_ARRAY_TASK_ID / 7))
q_index=$((SLURM_ARRAY_TASK_ID % 7))

# Extract parameters
s=${ss[$s_index]}
q=${qs[$q_index]}

echo "Running animator with s=$s and q=$q"

python "./Unity Analysis/Spring 2026/frac_map_animator.py" --mass_ratio $q --separation $s