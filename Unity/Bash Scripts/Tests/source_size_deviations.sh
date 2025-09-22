#!/usr/bin/env bash

#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=source_size_deviations
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --output="./Unity/Output Logs/Tests/source_size_deviations.txt"

# Commands to run
module load mamba
mamba activate .venv
python "./Unity/Python Scripts/Tests/source_size_deviations.py"