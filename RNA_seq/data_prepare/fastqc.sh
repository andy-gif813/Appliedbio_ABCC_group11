#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --output=fastqc_%j.log

module load fastqc/0.12.1-gcc-13.2.0

INDIR=/scratch/grp/msc_appbio/Group11_ABCC/RNA_seq/trimmed/R2
OUTDIR=/scratch/grp/msc_appbio/Group11_ABCC/RNA_seq/trimmed_fastqc
# mkdir -p "$OUTDIR"

for f in "$INDIR"/*.fastq; do
   echo "Processing $f"
   fastqc -o "$OUTDIR" "$f"
done
