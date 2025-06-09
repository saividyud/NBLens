#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --job-name=single_lens_1e11_cusp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu
#SBATCH --output="./Unity/Analysis 6-4/single_bash_logging2.txt"

# Commands to run
module load mamba
mamba activate .venv
python "./Unity/Analysis 6-4/single_lens_calculator2.py"

git add -A
git commit -m "Single lens 1e11 at cusp"
git push
