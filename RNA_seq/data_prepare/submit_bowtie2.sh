#!/bin/bash
# submit_bowtie2.sh - Submit bowtie2 alignment job

echo "=== Submitting RNA-seq alignment job ==="
echo "Current directory: $(pwd)"
echo "Time: $(date)"

# Check for active conda environment
if [[ -n "$CONDA_DEFAULT_ENV" ]]; then
    echo "Current conda environment: $CONDA_DEFAULT_ENV"
    # Deactivate conda environment to avoid conflicts with cluster module loading
    conda deactivate
fi

# Check if trimmed files exist
echo ""
echo "=== Checking input files ==="
TRIM_DIR="/scratch/grp/msc_appbio/Group11_ABCC/RNA_seq/trimmed"
if [ -d "$TRIM_DIR" ]; then
    echo "Found trimmed directory: $TRIM_DIR"
    echo "File list:"
    ls -lh "$TRIM_DIR"/*.fastq | head -12
else
    echo "Error: Cannot find trimmed directory $TRIM_DIR"
    exit 1
fi

# Check reference genome index
echo ""
echo "=== Checking reference genome index ==="
REF_DIR="/scratch/grp/msc_appbio/Group11_ABCC/RNA_seq/S288c"
if [ -d "$REF_DIR" ]; then
    echo "Reference genome directory: $REF_DIR"
    echo "Index files:"
    ls -lh "$REF_DIR"/*.bt2 2>/dev/null || echo "No .bt2 index files found"
    ls -lh "$REF_DIR"/*.bt2l 2>/dev/null || echo "No .bt2l index files found"
else
    echo "Warning: Reference genome directory not found $REF_DIR"
    echo "Please ensure correct index path is set in bowtie2 script"
fi

# Create working directory
WORK_DIR="/scratch/grp/msc_appbio/Group11_ABCC/RNA_seq/alignment_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

echo ""
echo "=== Creating working directory ==="
echo "Working directory: $WORK_DIR"

# Copy bowtie2 script to working directory
SCRIPT_DIR=$(dirname "$0")
BOWTIE2_SCRIPT="bowtie2_all.sh"

if [ -f "$BOWTIE2_SCRIPT" ]; then
    cp "$BOWTIE2_SCRIPT" "$WORK_DIR/"
    echo "Copied script: $BOWTIE2_SCRIPT"
else
    echo "Error: Cannot find bowtie2 script $BOWTIE2_SCRIPT"
    echo "Please ensure script is in current directory or specify full path"
    exit 1
fi

# Modify path in script (if needed)
sed -i "s|REFERENCE_INDEX=\".*\"|REFERENCE_INDEX=\"$REF_DIR/S288C_index\"|g" "$WORK_DIR/$BOWTIE2_SCRIPT"

# Submit job
echo ""
echo "=== Submitting Slurm job ==="
cd "$WORK_DIR"
sbatch "$BOWTIE2_SCRIPT"

echo ""
echo "=== Job submitted ==="
echo "Working directory: $WORK_DIR"
echo "Check job status: squeue -u $USER"
echo "Check output logs: $WORK_DIR/bowtie2_all_*.out"
