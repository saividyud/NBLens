#!/usr/bin/env bash

#SBATCH --time=44:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=single_lens_array_0.7
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-3

#SBATCH --output="./Unity/Output Logs/Single_Collection/single_0.7_%a_output.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
ss=0.7
qs=(1e-6 3e-6 1e-5 3e-5)

# Compute indices
q_index=$((SLURM_ARRAY_TASK_ID % 4))

# Extract parameters
s=$ss
q=${qs[$q_index]}

echo "Running $s separation with $q mass ratio with origin of (0, 0)"

python "./Unity/Python Scripts/Single Sims/single_binoffset.py" -s $s -q $q
