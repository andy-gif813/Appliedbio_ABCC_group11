#!/bin/bash
# Script: generate_new_slurm_scripts.sh
# Purpose: Batch generate Slurm submission scripts for RSEM analysis using new reference genome and trimmed data

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================

# Set working directory
WORKDIR="/scratch/grp/msc_appbio/Group11_ABCC/RNA_seq/Calculate_expression_level"

# Paths for RSEM analysis
NEW_INDEX="${WORKDIR}/S288c_R64/S288c_R64_index"     # New reference genome index path
RSEM_CMD="${WORKDIR}/RSEM/rsem-calculate-expression" # Local RSEM executable path
OUTPUT_DIR="${WORKDIR}/rsem_results_new"             # New output directory

# Directory containing trimmed FASTQ files
TRIMMED_DIR="../trimmed"

# ============================================================================
# SAMPLE MAPPING DEFINITION
# ============================================================================

# Define sample name to FASTQ prefix mapping (with your renamed BY4742_R1 files)
declare -A SAMPLE_MAP=(
    ["BY4742_R1"]="BY4742_R1_SRR7059707"
    ["BY4742_R2"]="BY4742_R2_SRR7059706"
    ["BY4742_R3"]="BY4742_R3_SRR7059709"
    ["SY14_R1"]="SY14_R1_SRR7059708"
    ["SY14_R2"]="SY14_R2_SRR7059705"
    ["SY14_R3"]="SY14_R3_SRR7059704"
)

# ============================================================================
# INITIALIZATION AND VALIDATION
# ============================================================================

echo "=== Generating Slurm Scripts for RSEM Analysis ==="
echo "Work directory: $WORKDIR"
echo "Reference index: $NEW_INDEX"
echo "RSEM command: $RSEM_CMD"
echo "Output directory: $OUTPUT_DIR"
echo "Trimmed data directory: $TRIMMED_DIR"
echo ""

# Change to working directory
cd "$WORKDIR" || { echo "Error: Cannot change to directory $WORKDIR"; exit 1; }

# Create necessary directories
echo "Creating necessary directories..."
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# Validate reference genome index
echo "Checking reference genome index..."
if [ ! -f "${NEW_INDEX}.1.bt2" ] && [ ! -f "${NEW_INDEX}.1.bt2l" ]; then
    echo "WARNING: Reference genome index not found or incomplete!"
    echo "Expected files: ${NEW_INDEX}.1.bt2 (or .bt2l) and other index files"
    echo "Please ensure the index is properly built."
else
    echo "✓ Reference genome index found"
fi

# Validate RSEM executable
echo "Checking RSEM executable..."
if [ ! -f "$RSEM_CMD" ]; then
    echo "ERROR: RSEM executable not found at: $RSEM_CMD"
    echo "Please check the path or install RSEM."
    exit 1
else
    echo "✓ RSEM executable found"
fi

# Validate trimmed directory
echo "Checking trimmed data directory..."
if [ ! -d "$TRIMMED_DIR" ]; then
    echo "ERROR: Trimmed data directory not found: $TRIMMED_DIR"
    exit 1
else
    echo "✓ Trimmed data directory found"
fi

# ============================================================================
# SCRIPT GENERATION
# ============================================================================

echo ""
echo "Generating Slurm scripts for the following samples:"

# Counter for successful script generation
GENERATED_COUNT=0
MISSING_FILES_COUNT=0

for SAMPLE in "${!SAMPLE_MAP[@]}"; do
    PREFIX="${SAMPLE_MAP[$SAMPLE]}"
    R1_FILE="${TRIMMED_DIR}/${PREFIX}_1_paired.fastq"
    R2_FILE="${TRIMMED_DIR}/${PREFIX}_2_paired.fastq"
    SCRIPT_NAME="submit_rsem_${SAMPLE}_new_ref.slurm"
    
    echo ""
    echo "Processing sample: $SAMPLE"
    echo "  FastQ prefix: $PREFIX"
    echo "  R1 file: $R1_FILE"
    echo "  R2 file: $R2_FILE"
    
    # Check if input files exist
    if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
        echo "  WARNING: Input file(s) missing! Skipping this sample."
        MISSING_FILES_COUNT=$((MISSING_FILES_COUNT + 1))
        continue
    fi
    
    # Get file sizes for logging
    R1_SIZE=$(ls -lh "$R1_FILE" 2>/dev/null | awk '{print $5}' || echo "unknown")
    R2_SIZE=$(ls -lh "$R2_FILE" 2>/dev/null | awk '{print $5}' || echo "unknown")
    
    echo "  File sizes: R1=$R1_SIZE, R2=$R2_SIZE"
    
    # Create the Slurm script
    cat > "$SCRIPT_NAME" << SCRIPT
#!/bin/bash
#SBATCH --job-name=rsem_${SAMPLE}
#SBATCH --output=${WORKDIR}/logs/rsem_${SAMPLE}_new_%j.out
#SBATCH --error=${WORKDIR}/logs/rsem_${SAMPLE}_new_%j.err
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=msc_appbio
#SBATCH --account=kcl

# ============================================================================
# JOB EXECUTION SECTION
# ============================================================================

echo "=== RSEM Analysis for Sample: ${SAMPLE} ==="
echo "Start time: \$(date)"
echo "Job ID: \$SLURM_JOB_ID"
echo "Sample: ${SAMPLE}"
echo "Input files:"
echo "  R1: ${R1_FILE}"
echo "  R2: ${R2_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Change to working directory
cd "${WORKDIR}" || { echo "Error: Cannot change to directory ${WORKDIR}"; exit 1; }

# Load required module
module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6

echo "Running RSEM analysis..."
echo "Command: ${RSEM_CMD} --bowtie2 --paired-end --estimate-rspd --append-names --num-threads 8 \\"
echo "  ${R1_FILE} \\"
echo "  ${R2_FILE} \\"
echo "  ${NEW_INDEX} \\"
echo "  ${OUTPUT_DIR}/${SAMPLE}"
echo ""

# Execute RSEM
${RSEM_CMD} \\
    --bowtie2 \\
    --paired-end \\
    --estimate-rspd \\
    --append-names \\
    --num-threads 8 \\
    "${R1_FILE}" \\
    "${R2_FILE}" \\
    "${NEW_INDEX}" \\
    "${OUTPUT_DIR}/${SAMPLE}"

# Capture exit status
EXIT_CODE=\$?
echo ""
echo "RSEM execution completed"
echo "Exit code: \$EXIT_CODE"
echo "End time: \$(date)"

# Check results
if [ \$EXIT_CODE -eq 0 ]; then
    echo "✓ RSEM analysis successful for ${SAMPLE}"
    
    # List generated files
    echo "Generated files:"
    ls -lh "${OUTPUT_DIR}/${SAMPLE}"* 2>/dev/null | head -5
    
    # Check gene count
    RESULT_FILE="${OUTPUT_DIR}/${SAMPLE}.genes.results"
    if [ -f "\$RESULT_FILE" ]; then
        LINES=\$(wc -l < "\$RESULT_FILE")
        GENES=\$((LINES-1))
        echo "Number of genes analyzed: \$GENES"
    fi
else
    echo "✗ RSEM analysis failed for ${SAMPLE}"
    echo "Please check the error log for details."
fi

echo "=== Job completed ==="
SCRIPT
    
    # Make the script executable
    chmod +x "$SCRIPT_NAME"
    
    echo "  ✓ Generated: $SCRIPT_NAME"
    GENERATED_COUNT=$((GENERATED_COUNT + 1))
done

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "=== Script Generation Summary ==="
echo "Total samples: ${#SAMPLE_MAP[@]}"
echo "Successfully generated: $GENERATED_COUNT"
echo "Skipped (missing files): $MISSING_FILES_COUNT"

if [ $GENERATED_COUNT -gt 0 ]; then
    echo ""
    echo "=== Generated Scripts ==="
    ls -l submit_rsem_*_new_ref.slurm 2>/dev/null
    
    echo ""
    echo "=== To Submit All Jobs ==="
    echo "You can submit all generated scripts using:"
    echo ""
    echo "  # Method 1: Submit individually"
    echo "  for script in submit_rsem_*_new_ref.slurm; do"
    echo "    echo \"Submitting \$script...\""
    echo "    sbatch \$script"
    echo "    sleep 2"
    echo "  done"
    echo ""
    echo "  # Method 2: Using a submission script"
    echo "  cat > submit_all_new.sh << 'EOF'"
    echo "  #!/bin/bash"
    echo "  echo \"Submitting all RSEM jobs...\""
    echo "  for script in submit_rsem_*_new_ref.slurm; do"
    echo "    echo \"Submitting: \$script\""
    echo "    JOB_ID=\$(sbatch --parsable \$script)"
    echo "    echo \"  Job ID: \$JOB_ID\""
    echo "    sleep 2"
    echo "  done"
    echo "  echo \"All jobs submitted!\""
    echo "  EOF"
    echo "  chmod +x submit_all_new.sh"
    echo "  ./submit_all_new.sh"
    echo ""
    echo "=== Important Notes ==="
    echo "1. Check that the reference genome index is built: ls -l ${NEW_INDEX}*"
    echo "2. Ensure trimmed data files exist in: $TRIMMED_DIR"
    echo "3. Output will be saved to: $OUTPUT_DIR"
    echo "4. Logs will be saved to: $WORKDIR/logs/"
else
    echo ""
    echo "ERROR: No scripts were generated!"
    echo "Please check:"
    echo "1. Trimmed data files exist in $TRIMMED_DIR"
    echo "2. File naming matches the pattern: {PREFIX}_1_paired.fastq and {PREFIX}_2_paired.fastq"
    echo "3. The SAMPLE_MAP array contains correct prefixes"
fi

echo ""
echo "=== Script generation completed ==="
echo "Time: $(date)"
