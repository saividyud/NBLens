#!/usr/bin/env bash

#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=single_lens_array
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-62

#SBATCH --output="./Unity/Output Logs/Single_Collection/single_%a_output.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
s1s=(0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
q1s=(1e-6 3e-6 1e-5 3e-5 1e-4 3e-4 1e-3)

# Compute indices
s1_index=$((SLURM_ARRAY_TASK_ID / 7))
q1_index=$((SLURM_ARRAY_TASK_ID % 7))

# Extract parameters
s1=${s1s[$s1_index]}
q1=${q1s[$q1_index]}

echo "Running $s1 separation with $q1 mass ratio with origin of binary offset"

python "./Unity/Python Scripts/Binary Sims/binary_binoffset.py" -s1 $s1 -q1 $q1
