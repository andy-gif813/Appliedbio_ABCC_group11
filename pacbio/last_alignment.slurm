cat > 02_scripts/last_alignment.slurm << 'EOF'
#!/bin/bash
#SBATCH --job-name=last_align
#SBATCH --output=03_logs/last_align_%j.log
#SBATCH --error=03_logs/last_align_%j.err
#SBATCH --time=12:00:00        # Reserve 12 hours, sufficient but not excessive queue usage
#SBATCH --mem=32G              # 32GB memory, enough for LAST
#SBATCH --cpus-per-task=8      # 8 CPU cores (even if LAST is not compiled with multithreading, other steps can use them)
#SBATCH --nodes=1

echo "======= LAST Whole Genome Alignment Task ======="
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo "====================================="

# 1. Robust environment activation (using previously verified successful method)
source ~/.bashrc  # Load conda initialization
conda activate genome_comparison

echo "=== Environment Verification ==="
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "LAST version: $(lastal --version 2>&1 | head -1)"
echo "Python version: $(python --version 2>&1)"

# 2. Define file paths
SY14_ASSEMBLY="04_results/flye_assembly_full/assembly.fasta"
REF_GENOME="01_reference_genome/yeast_genome.fna"
OUT_DIR="04_results/genome_comparison"
OUT_PREFIX="${OUT_DIR}/sy14_vs_by4742"

mkdir -p "$OUT_DIR"

# 3. Check input files
echo ""
echo "=== Input File Verification ==="
if [ ! -f "$SY14_ASSEMBLY" ]; then
    echo "❌ Error: SY14 assembly file does not exist: $SY14_ASSEMBLY"
    exit 1
fi
if [ ! -f "$REF_GENOME" ]; then
    echo "❌ Error: Reference genome file does not exist: $REF_GENOME"
    exit 1
fi
echo "✅ SY14 assembly: $(ls -lh $SY14_ASSEMBLY)"
echo "✅ Reference genome: $(ls -lh $REF_GENOME)"

# 4. Build LAST database (one of the most time-consuming steps)
echo ""
echo "=== Step 1/3: Building LAST Database ==="
echo "Building reference genome database..."
SECONDS=0
lastdb -v "${OUT_PREFIX}_ref_db" "$REF_GENOME"
echo "Reference database building completed, time taken: $((SECONDS/60)) minutes $((SECONDS%60)) seconds"

echo "Building query sequence database..."
SECONDS=0
lastdb -v "${OUT_PREFIX}_qry_db" "$SY14_ASSEMBLY"
echo "Query database building completed, time taken: $((SECONDS/60)) minutes $((SECONDS%60)) seconds"

# 5. Run LAST alignment
echo ""
echo "=== Step 2/3: Running LAST Whole Genome Alignment ==="
echo "Starting alignment..."
SECONDS=0

# Use more appropriate parameters for yeast genome alignment
lastal \
    -v \               # Show progress information
    -e 0.05 \          # E-value threshold
    -a 1 \             # Match score
    -b 1 \             # Mismatch penalty
    "${OUT_PREFIX}_ref_db" "$SY14_ASSEMBLY" > "${OUT_PREFIX}.maf"

ALIGN_SECONDS=$SECONDS
echo "Alignment completed, total time taken: $((ALIGN_SECONDS/60)) minutes $((ALIGN_SECONDS%60)) seconds"

# 6. Result verification and statistics
echo ""
echo "=== Step 3/3: Result Verification ==="
if [ -f "${OUT_PREFIX}.maf" ] && [ -s "${OUT_PREFIX}.maf" ]; then
    echo "✅ Alignment successful!"
    echo "   Output file: ${OUT_PREFIX}.maf"
    echo "   File size: $(ls -lh ${OUT_PREFIX}.maf | awk '{print $5}')"
    echo "   Number of lines: $(wc -l < ${OUT_PREFIX}.maf)"
    
    # Extract brief statistical information
    echo ""
    echo "=== Alignment Result Preview ==="
    echo "First 5 alignment blocks:"
    grep -c "^a score" "${OUT_PREFIX}.maf" | awk '{print "   Total alignment blocks: "$1}'
    echo ""
    echo "MAF file header example:"
    head -20 "${OUT_PREFIX}.maf"
else
    echo "❌ Error: No valid MAF file generated"
    exit 1
fi

# 7. Automatically convert to Circos format (optional but recommended)
echo ""
echo "=== Extra Step: Converting to Circos Format ==="
if [ -f "02_scripts/maf_to_circos.py" ]; then
    echo "Conversion script found, generating Circos link file..."
    python3 02_scripts/maf_to_circos.py \
        "${OUT_PREFIX}.maf" \
        "${OUT_DIR}/circos_links.txt"
    
    if [ -f "${OUT_DIR}/circos_links.txt" ]; then
        echo "✅ Circos link file generated: ${OUT_DIR}/circos_links.txt"
        echo "   Number of links: $(wc -l < ${OUT_DIR}/circos_links.txt)"
        echo "   First 3 links:"
        head -3 "${OUT_DIR}/circos_links.txt"
    fi
else
    echo "⚠️  maf_to_circos.py script not found, please convert manually."
fi

echo ""
echo "======= Task Completed ======="
echo "End time: $(date)"
echo "Output directory: $OUT_DIR"
echo "Main files:"
ls -lh "$OUT_DIR"/*.maf "$OUT_DIR"/*.txt 2>/dev/null | head -10
EOF

# Submit the job
sbatch 02_scripts/last_alignment.slurm
