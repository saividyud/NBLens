#!/usr/bin/env bash
# CPU array job: runs 125 simulations with up to 40 in parallel.
# Complements triple_binoffset_gpu.sh; both skip already-completed sims
# via lock files so no work is duplicated.
#
# Submit from project root: sbatch "Unity/Bash Scripts/Part 5/triple_binoffset_cpu.sh"

#SBATCH --time=25:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=triple_cpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-124

#SBATCH --output="./Unity/Output Logs/Part_5/triple_cpu_%a_2.txt"

# Commands to run
module load mamba
mamba activate .venv

python "./Unity/Python Scripts/Part 5/triple_binoffset.py" \
    --task-id $SLURM_ARRAY_TASK_ID --mode cpu
