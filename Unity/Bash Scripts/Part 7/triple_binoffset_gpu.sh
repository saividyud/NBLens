#!/usr/bin/env bash
# GPU series job: churns through all 60 sims in REVERSE order (59 -> 0).
# Complements triple_binoffset_cpu.sh; both skip already-completed sims
# via lock files so no work is duplicated.
#
# Submit from project root: sbatch "Unity/Bash Scripts/Part 7/triple_binoffset_gpu.sh"

#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=triple_gpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --gres=gpu:1

#SBATCH --output="./Unity/Output Logs/Part_7/triple_gpu_%j.txt"

# Commands to run
module load mamba
module load cuda/12.9.0-vvrwkrx
export LD_LIBRARY_PATH=$CUDA_HOME/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
mamba activate .venv

rm -rf ~/.cupy/kernel_cache/

python "./Unity/Python Scripts/Part 7/triple_binoffset.py" \
    --batch-start 0 --batch-end 60 --mode gpu --reverse
