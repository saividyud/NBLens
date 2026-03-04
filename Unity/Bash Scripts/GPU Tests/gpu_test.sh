#!/usr/bin/env bash
# GPU allocation test script for Unity cluster
# Tests that SLURM assigns a GPU node and CuPy can access it.
#
# Submit from project root: sbatch "Unity/Bash Scripts/GPU Tests/gpu_test.sh"

#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=gpu_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

# Request 1 GPU (Unity has 16 NVIDIA GPUs across ~107 nodes)
#SBATCH --gres=gpu:1

#SBATCH --output="./Unity/Output Logs/GPU_Tests/gpu_test_%j_output.txt"

# Commands to run
module load mamba
mamba load cuda/12.9.0-vvrwkrx
mamba activate .venv

# Ensure output directory exists
mkdir -p "./Unity/Output Logs/GPU_Tests"

echo "=== SLURM Job Info ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Allocated GPUs: $CUDA_VISIBLE_DEVICES"
echo ""

python "./Unity/Python Scripts/GPU Tests/gpu_node_stats.py"
