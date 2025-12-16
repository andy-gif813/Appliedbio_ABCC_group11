#!/bin/bash
#SBATCH --job-name=indel_calling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=indel_calling_%j.log

module load samtools/1.17
module load bcftools/1.12

# create mpileup
bcftools mpileup -f WT_genome.fasta SY14.sorted.bam -Ou -a FORMAT/DP > SY14.mpileup.bcf

# call variants (SNP + indel)
bcftools call -mv -Oz -o SY14.raw.vcf.gz SY14.mpileup.bcf

# VCF
bcftools index SY14.raw.vcf.gz

# filter other things and left indel
bcftools view -v indels SY14.raw.vcf.gz -Oz -o SY14.indel.vcf.gz
bcftools index SY14.indel.vcf.gz

echo "Indel calling completed!"
