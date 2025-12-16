#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --output=trim_%j.log

# Load Trimmomatic module
module load trimmomatic

# adaptor file
ADAPTER_FILE="/software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/trimmomatic-0.39-i5lbieeefpcxmue5vf4sfmh2p636cj27/share/adapters/TruSeq3-PE.fa"

# Input parent directory
INDIR=/scratch/grp/msc_appbio/Group11_ABCC/RNA_seq/rawdata/BY4742_R2

# Output directory
OUTDIR=/scratch/grp/msc_appbio/Group11_ABCC/RNA_seq/trimmed

# Loop and Run Trimmomatic (PE mode)
for R1 in $INDIR/*_1.fastq; do
	SAMPLE=$(basename "$R1" | sed 's/_1.fastq//')
	R2="$INDIR/${SAMPLE}_2.fastq"

	echo "Processing sample: $SAMPLE"

	## Output file names
    OUT_R1_PAIRED="$OUTDIR/${SAMPLE}_1_paired.fastq"
    OUT_R1_UNPAIRED="$OUTDIR/${SAMPLE}_1_unpaired.fastq"
    OUT_R2_PAIRED="$OUTDIR/${SAMPLE}_2_paired.fastq"
    OUT_R2_UNPAIRED="$OUTDIR/${SAMPLE}_2_unpaired.fastq"


	trimmomatic PE -threads 8 \
    "$R1" "$R2" \
    "$OUT_R1_PAIRED" \
    "$OUT_R1_UNPAIRED" \
    "$OUT_R2_PAIRED" \
    "$OUT_R2_UNPAIRED"  \
    ILLUMINACLIP:$ADAPTER_FILE:2:30:10 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    HEADCROP:16
	
   echo "Finished trimming $SAMPLE"
done

