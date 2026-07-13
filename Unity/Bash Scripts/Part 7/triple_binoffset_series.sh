#!/usr/bin/env bash
# Batched triple lens simulations on GPU.
# 60 total simulations split across SIMS_PER_JOB simulations per job.
# Submit from project root: sbatch "Unity/Bash Scripts/Part 7/triple_binoffset_series.sh"

# ── Tunable: how many simulations each GPU job runs ──
SIMS_PER_JOB=5

#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=triple_lens_series
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

# 60 sims / 5 per job = 12 jobs (array indices 0-11)
#SBATCH --array=0-11
#SBATCH --gres=gpu:1

#SBATCH --output="./Unity/Output Logs/Part_7/triple_series_%a.txt"

# Commands to run
module load mamba
module load cuda/12.9.0-vvrwkrx
export LD_LIBRARY_PATH=$CUDA_HOME/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
mamba activate .venv

BATCH_START=$((SLURM_ARRAY_TASK_ID * SIMS_PER_JOB))
BATCH_END=$((BATCH_START + SIMS_PER_JOB))

echo "Job $SLURM_ARRAY_TASK_ID: running simulations $BATCH_START to $((BATCH_END - 1))"

python "./Unity/Python Scripts/Part 7/triple_binoffset_series.py" \
    --batch-start $BATCH_START \
    --batch-end $BATCH_END
