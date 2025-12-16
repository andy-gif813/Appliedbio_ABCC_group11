#!/bin/bash
#SBATCH --job-name=indel_filter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=indel_filter_%j.log

module load bcftools/1.12
module load htslib/1.12

# Step 1: filter indel
# QUAL>=30, depth FORMAT/DP>=10 and FORMAT/DP<=500, IMF>=0.3
bcftools filter -i 'QUAL>=30 && FORMAT/DP>=10 && FORMAT/DP<=500 && IMF>=0.3 && TYPE="indel"' \
    SY14.indel.vcf.gz -Oz -o SY14.indel.filtered.vcf.gz

# Step 2: creat
bcftools index SY14.indel.filtered.vcf.gz

# Step 3: Count Indel
# Count indel only, not header
INDL_COUNT=$(bcftools view -i 'TYPE="indel"' SY14.indel.filtered.vcf.gz | grep -v "^#" | wc -l)

echo "Filtered indel count: $INDL_COUNT"

# === Done ===
echo "Indel filtering and counting completed!"
