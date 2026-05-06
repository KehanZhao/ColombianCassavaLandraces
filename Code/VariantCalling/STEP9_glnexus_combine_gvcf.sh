#!/bin/bash -l
#SBATCH -J glnexus
#SBATCH -t 148:00:00
#SBATCH -c 48
#SBATCH --mem=1048G
#SBATCH --ntasks=1
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/glnexus-%j.out

cd ./VCF_withKistler_noHerbarium

# Module load
module load apptainer

# Running apptainer to build the docker environment
apptainer run docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1 \
glnexus_cli \
-c DeepVariant \
../deepvariant_output/*.gvcf.gz > combined_cassava_mar2025_final_withKistler_noHerbarium.bcf

echo "glnexus completed"

# Convert bcf into gvcf

module load bcftools htslib

bcftools view combined_cassava_mar2025_final_withKistler_noHerbarium.bcf | bgzip -@ 16 -c > combined_cassava_mar2025_final_withKistler_noHerbarium.vcf.gz
