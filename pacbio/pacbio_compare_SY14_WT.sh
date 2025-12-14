#!/bin/bash
#SBATCH --job-name=correct_path_analysis
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=03_logs/correct_path_analysis_%j.log

echo "=== Analysis Using Correct Paths ==="
echo "Start time: $(date)"
echo ""

# Load system modules
module purge
module load samtools/1.17-gcc-13.2.0-python-3.11.6 2>/dev/null

# Use minimap2 from conda environment
MINIMAP2_PATH="/scratch/grp/msc_appbio/Group11_ABCC/PacBio/05_envs/flye_env/bin/minimap2"

echo "1. Tool verification:"
echo "minimap2: $($MINIMAP2_PATH --version 2>&1 | head -1)"
echo "samtools: $(samtools --version 2>&1 | head -1)"

# Set working directory to project root
cd /scratch/grp/msc_appbio/Group11_ABCC/PacBio
echo -e "\nProject root directory: $(pwd)"

# Create output directory
OUTPUT="04_results/correct_path_analysis_$(date +%Y%m%d_%H%M%S)"
mkdir -p $OUTPUT
cd $OUTPUT

echo "Output directory: $(pwd)"
echo ""

# Define file paths (relative to project root)
REF="/scratch/grp/msc_appbio/Group11_ABCC/PacBio/01_reference_genome/yeast_genome.fna"
DATA="/scratch/grp/msc_appbio/Group11_ABCC/PacBio/00_raw_data/SRR6823435_subreads.fastq.gz"

echo "Reference genome: $REF"
echo "Data file: $DATA"
echo ""

# 2. Extract small data subset (0.1%)
echo "2. Extracting data subset (0.1%)..."
SAMPLE_READS=500
echo "Sampling: $SAMPLE_READS reads"

# Check if files exist
if [ ! -f "$DATA" ]; then
  echo "Error: Data file not found: $DATA"
  exit 1
fi

if [ ! -f "$REF" ]; then
  echo "Error: Reference genome not found: $REF"
  exit 1
fi

zcat $DATA | head -$((SAMPLE_READS * 4)) > sy14_sample.fastq

echo "Sample statistics:"
READ_COUNT=$(grep -c '^@' sy14_sample.fastq)
echo "Number of reads: $READ_COUNT"

if [ $READ_COUNT -eq 0 ]; then
  echo "Error: No reads extracted"
  exit 1
fi

AVG_LENGTH=$(awk 'NR%4==2{sum+=length($1); count++} END{if(count>0)print sum/count; else print 0}' sy14_sample.fastq)
echo "Average read length: $AVG_LENGTH bp"
echo "File size: $(ls -lh sy14_sample.fastq | awk '{print $5}')"
echo ""

# 3. Align to reference genome
echo "3. Aligning SY14 reads to BY4742 reference genome..."
$MINIMAP2_PATH -ax map-pb -t 16 $REF sy14_sample.fastq > alignment.sam 2> minimap2.log

echo "Alignment completed, checking results..."
if [ -s "alignment.sam" ]; then
  echo "✓ SAM file created successfully"
  echo "Number of lines in SAM file: $(wc -l < alignment.sam)"
  echo "File size: $(ls -lh alignment.sam | awk '{print $5}')"
else
  echo "✗ SAM file is empty, check log:"
  cat minimap2.log
  exit 1
fi

echo ""

# 4. Process alignment results
echo "4. Processing alignment results..."
samtools view -bS alignment.sam > alignment.bam
samtools sort alignment.bam -o sorted.bam
samtools index sorted.bam

echo "5. Alignment statistics:"
samtools flagstat sorted.bam > flagstat.txt
cat flagstat.txt

echo -e "\n6. Coverage analysis:"
samtools depth sorted.bam > coverage.txt

if [ -s "coverage.txt" ]; then
  # Reference genome size
  REF_SIZE=$(grep -v '^>' $REF | tr -d '\n' | wc -c)
  TOTAL_POSITIONS=$(wc -l < coverage.txt)
  COV_POSITIONS=$(awk '$3>0' coverage.txt | wc -l)
  AVG_DEPTH=$(awk '{sum+=$3} END{print sum/NR}' coverage.txt)

  echo "Reference genome size: $REF_SIZE bp"
  echo "Number of analyzed positions: $TOTAL_POSITIONS"
  echo "Number of covered positions: $COV_POSITIONS"
  echo "Average depth: $AVG_DEPTH"

  if [ $REF_SIZE -gt 0 ]; then
    COVERAGE_PERCENT=$(echo "scale=2; $COV_POSITIONS * 100 / $REF_SIZE" | bc 2>/dev/null || echo "Calculation failed")
    echo "Genome coverage: $COVERAGE_PERCENT%"
  fi
else
  echo "Coverage file is empty"
fi

echo -e "\n7. Generating result summary..."
cat > analysis_summary.txt << SUMMARY
SY14 vs BY4742 Comparison Analysis
=======================
Analysis time: $(date)
Data sample: $READ_COUNT reads
Reference genome: BY4742

Alignment statistics:
$(cat flagstat.txt)

Coverage statistics:
- Number of analyzed positions: $TOTAL_POSITIONS
- Number of covered positions: $COV_POSITIONS
- Average depth: $AVG_DEPTH
- Genome coverage: ${COVERAGE_PERCENT}%

Analysis description:
SY14 PacBio long-read data was directly aligned to the BY4742 reference genome.
This analysis shows preliminary alignment between the two strains.

Next steps:
1. Increase data volume for better statistical significance
2. Perform deeper analysis with more reads
3. Identify specific variant positions

SUMMARY

echo -e "\n8. Organizing result files..."
mkdir -p results
mv alignment.sam alignment.bam sorted.bam sorted.bam.bai coverage.txt flagstat.txt minimap2.log results/ 2>/dev/null
cp analysis_summary.txt results/

echo -e "\nAnalysis completed!"
echo "Results saved in: $(pwd)"
echo "Key files:"
ls -lh results/
echo ""
echo "End time: $(date)"
