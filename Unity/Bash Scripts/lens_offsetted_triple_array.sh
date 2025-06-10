#!/usr/bin/env bash

#SBATCH --time=60:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=triple_lens_offsetted_array
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-14

#SBATCH --output="./Unity/Output Logs/Collection_0.8/triple_1e11_%a_output.txt"
#SBATCH --error="./Unity/Error Logs/Collection_0.8/triple_1e11_%a_error.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
alphas=(0 45 90 135 180)
pmrs=(0.1 0.01 0.001)

# Compute indices
A_index=$(( SLURM_ARRAY_TASK_ID / 4 ))
S_index=$(( SLURM_ARRAY_TASK_ID % 4 ))

echo $A_index
echo $S_index

