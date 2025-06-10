#!/usr/bin/env bash

#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=lens_collection_0.8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --output="./Unity/Output Logs/Collection_0.8/single_1e11_%x_%s_output.txt"
#SBATCH --error="./Unity/Error Logs/Collection_0.8/single_1e11_%x_%s_error.txt"

# Commands to run
module load mamba
mamba activate .venv
python "./Unity/Python Scripts/lenses_offsetted.py" -a2 45 -pmr 1e-3 -l single