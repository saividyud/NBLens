#!/usr/bin/env bash
# GPU benchmark timing script
# Submit from project root: sbatch "Unity Analysis/Part 5/gpu_benchmark.sh"

#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=gpu_benchmark
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --gres=gpu:1

#SBATCH --output="./Unity/Output Logs/Part_5/gpu_benchmark_output_%j.txt"

# Commands to run
module load mamba
module load cuda/12.9.0-vvrwkrx
export LD_LIBRARY_PATH=$CUDA_HOME/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
mamba activate .venv

rm -rf ~/.cupy/kernel_cache/

echo "Running GPU benchmark"

python "./Unity Analysis/Part 5/gpu_benchmark.py" --warmup 1 --trials 15
