#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --job-name=single_lens_1e11_cusp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --output="./Unity/Output Logs/single_bash2_output.txt"
#SBATCH --error="./Unity/Error Logs/single_bash2_error.txt"

# Commands to run
module load mamba
mamba activate .venv
python "./Unity/Python Scripts/single_lens_calculator2.py"

git add -A
git commit -m "Single lens 1e11 at cusp"
git push
