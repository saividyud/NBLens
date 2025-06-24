#!/usr/bin/env bash

#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=varried_s_triple_offsetted_array
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-32

#SBATCH --output="./Unity/Output Logs/Collection_pmr_0.01/s_varried_%a_output.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
seperations=(0.5 0.6 0.7 0.8 0.9 1 1/0.9 1/0.8 1/0.7 1/0.6 1/0.5)
offsets=(cm binary_offset triple_offset)

# Compute indices
seperations_index=$((SLURM_ARRAY_TASK_ID / 3))
offsets_index=$((SLURM_ARRAY_TASK_ID % 3))

# Extract parameters
seperation=${seperations[$seperations_index]}
offset=${offsets[$offsets_index]}

echo "Running $alpha degrees with $pmr mass ratio and $seperation seperation with $offset offset"

python "./Unity/Python Scripts/lenses_offsetted.py" -s2 $seperation -a2 45 -pmr 0.01 -l triple -o $offset
