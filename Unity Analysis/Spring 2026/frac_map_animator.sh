#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=binary_animator_3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-4

#SBATCH --output="./Unity/Output Logs/Analysis/animator_%a_output.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
rs=(1e-1 3e-1 1e0 3e0 1e1)

r=${rs[$SLURM_ARRAY_TASK_ID]}

echo "Running animator with radius multiplier of $r"

python "./Unity Analysis/Spring 2026/frac_map_animator.py" --multiplier $r