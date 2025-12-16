#!/bin/bash
#SBATCH --job-name=SY14Illunina_mapping_SY14assemble
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=Illunina_mapping_SY14assemble_%j.log


module load anaconda3/2022.10
source activate bioinfo
### Input files
R1="SRR6825081_1.fastq"
R2="SRR6825081_2.fastq"
REF="SY14_assemble.fasta"
THREADS=$SLURM_CPUS_PER_TASK

echo "Step 1: Index reference genome"
bwa index $REF
echo "Step 2: Mapping with bwa mem"
bwa mem -t $THREADS \
    -R "@RG\tID:SRR6825081\tSM:SRR6825081\tPL:ILLUMINA" \
    $REF \
    $R1 $R2 | samtools sort -@ $THREADS -o SRR6825081.sorted.bam

echo "Step 3: Index BAM"
samtools index SRR6825081.sorted.bam

echo "Step 4: Mapping QC"
samtools flagstat SRR6825081.sorted.bam > SRR6825081.flagstat.txt

echo "Step 5: Depth calculation"
samtools depth -a SRR6825081.sorted.bam > SRR6825081.depth.txt

echo "===== DONE! ====="
echo "Output files:"
echo " -  $R1 / $R2"
echo " - SRR6825081.sorted.bam"
echo " - SRR6825081.sorted.bam.bai"
echo " - SRR6825081.flagstat.txt"
echo " - SRR6825081.depth.txt"
