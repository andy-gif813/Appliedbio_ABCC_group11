#!/bin/bash
# Script: rerun_all_origin.sh
# Purpose: Re-run RSEM analysis for all 6 original data samples

WORKDIR="/scratch/grp/msc_appbio/Group11_ABCC/Calculate_expression_level"
cd "$WORKDIR"

# Clean output directory
echo "Cleaning output directory..."
rm -rf rsem_results_origin
mkdir -p rsem_results_origin
mkdir -p logs

# Create log file
LOG_FILE="logs/rerun_all_origin_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================="
echo "Starting re-run of RSEM analysis for all 6 samples"
echo "Time: $(date)"
echo "Working directory: $WORKDIR"
echo "Output directory: rsem_results_origin"
echo "=========================================="

# SRR to sample mapping
declare -A SAMPLE_SRR_MAP=(
    ["BY4742_R1"]="SRR7059707"
    ["BY4742_R2"]="SRR7059706"
    ["BY4742_R3"]="SRR7059709"
    ["SY14_R1"]="SRR7059708"
    ["SY14_R2"]="SRR7059705"
    ["SY14_R3"]="SRR7059704"
)

# Raw data directory
RAW_DATA_BASE="/scratch/grp/msc_appbio/Group11_ABCC/rawdata"

# Generate Slurm scripts for each sample
for SAMPLE in "${!SAMPLE_SRR_MAP[@]}"; do
    SRR_ID="${SAMPLE_SRR_MAP[$SAMPLE]}"
    
    cat > "run_${SAMPLE}.slurm" << SCRIPT
#!/bin/bash
#SBATCH --job-name=origin_${SAMPLE}
#SBATCH --output=logs/origin_${SAMPLE}_%j.out
#SBATCH --error=logs/origin_${SAMPLE}_%j.err
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=msc_appbio
#SBATCH --account=kcl

module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6

cd "$WORKDIR"

echo "Starting processing ${SAMPLE}: \$(date)"
echo "Sample: ${SAMPLE}"
echo "SRR ID: ${SRR_ID}"
echo "Input files:"
echo "  R1: ${RAW_DATA_BASE}/${SAMPLE}/${SRR_ID}_1.fastq"
echo "  R2: ${RAW_DATA_BASE}/${SAMPLE}/${SRR_ID}_2.fastq"

# Check if input files exist
if [ ! -f "${RAW_DATA_BASE}/${SAMPLE}/${SRR_ID}_1.fastq" ]; then
    echo "Error: R1 file not found!"
    exit 1
fi

if [ ! -f "${RAW_DATA_BASE}/${SAMPLE}/${SRR_ID}_2.fastq" ]; then
    echo "Error: R2 file not found!"
    exit 1
fi

./RSEM/rsem-calculate-expression \\
    --bowtie2 \\
    --paired-end \\
    --estimate-rspd \\
    --append-names \\
    --num-threads 8 \\
    "${RAW_DATA_BASE}/${SAMPLE}/${SRR_ID}_1.fastq" \\
    "${RAW_DATA_BASE}/${SAMPLE}/${SRR_ID}_2.fastq" \\
    S288c_R64/S288c_R64_index \\
    rsem_results_origin/${SAMPLE}

echo "Completed ${SAMPLE}: \$(date)"
SCRIPT

    echo "Generated script: run_${SAMPLE}.slurm"
done

# Submit all jobs
echo ""
echo "Submitting all jobs..."
JOB_IDS=""

for SAMPLE in "${!SAMPLE_SRR_MAP[@]}"; do
    SCRIPT_FILE="run_${SAMPLE}.slurm"
    
    if [ -f "$SCRIPT_FILE" ]; then
        echo "Submitting: $SCRIPT_FILE"
        JOB_ID=$(sbatch --parsable "$SCRIPT_FILE")
        JOB_IDS="$JOB_IDS $JOB_ID"
        echo "  Job ID: $JOB_ID"
        sleep 2  # Avoid submitting too many jobs at once
    else
        echo "Warning: Script file not found: $SCRIPT_FILE"
    fi
done

echo "=========================================="
echo "All jobs submitted!"
echo "Job ID list:$JOB_IDS"
echo ""
echo "Use these commands to monitor job status:"
echo "  squeue -u \$USER"
echo "  watch -n 10 'squeue -u \$USER'"
echo ""
echo "View log files:"
echo "  ls -lt logs/origin_*.out"
echo ""
echo "Check result files:"
echo "  ls -la rsem_results_origin/"
echo "=========================================="
