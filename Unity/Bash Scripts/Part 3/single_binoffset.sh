#!/usr/bin/env bash

#SBATCH --time=44:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=single_lens_array_2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-1

#SBATCH --output="./Unity/Output Logs/Part_3/single_%a_output.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
qs=(3e-4 1e-3)

# Compute indices
q_index=$(SLURM_ARRAY_TASK_ID)

# Extract parameters
q=${qs[$q_index]}

echo "Running $s separation with $q mass ratio with origin of (0, 0)"

python "./Unity/Python Scripts/Part 3/single_binoffset.py" -s 0.9 -q $q
