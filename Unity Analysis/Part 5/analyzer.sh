#!/usr/bin/env bash

#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=triple_analyzer_5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-4

#SBATCH --output="./Unity/Output Logs/Analysis/analysis_%a_output.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
rs=(1e-1 3e-1 1e0 3e0 1e1)

r=${rs[$SLURM_ARRAY_TASK_ID]}

echo "Running analysis with radius multiplier of $r"

python "./Unity Analysis/Part 5/analyzer.py" --multiplier $r