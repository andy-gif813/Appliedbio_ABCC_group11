#!/bin/bash
#SBATCH --job-name=flye_full
#SBATCH --output=03_logs/flye_full_%j.log
#SBATCH --error=03_logs/flye_full_%j.err
#SBATCH --partition=cpu
#SBATCH --time=2-00:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1

# ====== 1. Robust Conda Environment Activation ======
# Method: Directly source conda configuration file, then activate environment (bypassing `conda activate`)
source /users/k25031516/miniconda3/etc/profile.d/conda.sh
conda activate flye_env

# Simple environment validation
echo "=== Environment Verification ==="
echo "Flye path: $(which flye)"
echo "Python path: $(which python)"

# ====== 2. Run Flye (first run, without --resume) ======
echo ""
echo "=== Starting Flye Whole Genome Assembly ==="

flye \
    --pacbio-raw 00_raw_data/SRR6823435_subreads.fastq.gz \
    --out-dir 04_results/flye_assembly_full \
    --threads 16
    # Note: Removed `--resume` parameter for the first run. Same as above.
