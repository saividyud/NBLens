#!/usr/bin/env bash

#SBATCH --time=44:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=binary_lens_array_2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-33

#SBATCH --output="./Unity/Output Logs/Part_3/binary_%a_output.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
ss=(0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1/0.9 1/0.8 1/0.7 1/0.6 1/0.5 1/0.4 1/0.3 1/0.2)
qs=(3e-4 1e-3)

# Compute indices
s_index=$((SLURM_ARRAY_TASK_ID / 2))
q_index=$((SLURM_ARRAY_TASK_ID % 2))

# Extract parameters
s=${ss[$s_index]}
q=${qs[$q_index]}

echo "Running $s separation with $q mass ratio with origin of binary offset"

python "./Unity/Python Scripts/Part 3/binary_binoffset.py" -s $s -q $q
