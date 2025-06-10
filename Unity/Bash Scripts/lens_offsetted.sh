#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --job-name=single_lens_1e11_offsetted_45deg
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --output="./Unity/Output Logs/single_lens_1e11_offsetted_45deg_output.txt"
#SBATCH --error="./Unity/Error Logs/single_lens_1e11_offsetted_45deg_error.txt"

# Commands to run
module load mamba
mamba activate .venv
python "./Unity/Python Scripts/single_lens_offsetted.py" -a2 45 -pmr 1e-3 -l single