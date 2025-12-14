#!/bin/bash
#SBATCH --job-name=bowtie2_all
#SBATCH --output=bowtie2_all_%j.out
#SBATCH --error=bowtie2_all_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu

echo "=== Batch Processing All Samples ==="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"

module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6
module load samtools/1.17-gcc-13.2.0

# Sample list
SAMPLES=("BY4742_R1" "BY4742_R2" "BY4742_R3" "SY14_R1" "SY14_R2" "SY14_R3")

# Create directories
mkdir -p aligned
mkdir -p logs

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "=== Processing $SAMPLE ==="
    
    # Find files
    R1=$(ls ../trimmed/${SAMPLE}_*_1_paired.fastq 2>/dev/null | head -1)
    R2=$(ls ../trimmed/${SAMPLE}_*_2_paired.fastq 2>/dev/null | head -1)
    
    if [ -z "$R1" ] || [ -z "$R2" ]; then
        echo "Error: Cannot find files for $SAMPLE"
        continue
    fi
    
    echo "Input files:"
    echo "  R1: $(basename $R1)"
    echo "  R2: $(basename $R2)"
    echo "  File size: R1=$(ls -lh $R1 | awk '{print $5}'), R2=$(ls -lh $R2 | aw>
    
    # Run bowtie2
    bowtie2 --threads 8 \
        -x S288c/S288C_index \
        -1 "$R1" \
        -2 "$R2" \
        -S aligned/${SAMPLE}.sam \
        2> logs/${SAMPLE}_bowtie2.log
    
    BOWTIE2_EXIT=$?
    if [ $BOWTIE2_EXIT -ne 0 ]; then
        echo "✗ bowtie2 failed, exit code: $BOWTIE2_EXIT"
        echo "Last 10 lines of error log:"
        tail -10 logs/${SAMPLE}_bowtie2.log
        continue
    fi
    
    echo "✓ bowtie2 successful"
    
    # Show alignment rate
    echo "Alignment rate:"
    grep "overall alignment rate" logs/${SAMPLE}_bowtie2.log
    
    # Convert and sort
    echo "Converting SAM to BAM..."
    samtools view -@ 15 -bS aligned/${SAMPLE}.sam > aligned/${SAMPLE}.bam
    
    echo "Sorting BAM..."
    samtools sort -@ 15 -o aligned/${SAMPLE}_sorted.bam aligned/${SAMPLE}.bam
    
    echo "Creating index..."
    samtools index -@ 15 aligned/${SAMPLE}_sorted.bam
    
    # Clean intermediate files
    rm aligned/${SAMPLE}.bam
    
    # Statistics information
    echo "BAM file information:"
    ls -lh aligned/${SAMPLE}_sorted.bam*
    
    echo "✓ Completed $SAMPLE"
done

echo ""
echo "=== Generating Summary Report ==="
echo "Sample,Total reads,Alignment rate,Properly paired rate" > logs/alignment_summary.csv
for SAMPLE in "${SAMPLES[@]}"; do
    if [ -f "aligned/${SAMPLE}_sorted.bam" ]; then
        # Get total reads count
        total=$(samtools view -c aligned/${SAMPLE}_sorted.bam)
        # Get alignment rate (from bowtie2 log)
        align_rate=$(grep "overall alignment rate" logs/${SAMPLE}_bowtie2.log |>
        # Get properly paired rate
        proper_paired=$(samtools flagstat aligned/${SAMPLE}_sorted.bam | grep ">
        proper_rate=$(echo "scale=2; $proper_paired * 100 / $total" | bc)
        
        echo "$SAMPLE,$total,$align_rate,${proper_rate}%" >> logs/alignment_sum>
    fi
done

echo ""
echo "=== Batch Processing Completed ==="
echo "End time: $(date)"
echo "Output directory: aligned/"
echo "Log directory: logs/"
echo "Summary report: logs/alignment_summary.csv"
