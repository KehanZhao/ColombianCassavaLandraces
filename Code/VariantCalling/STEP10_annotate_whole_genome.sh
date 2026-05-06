#!/bin/bash -l
#SBATCH -J rename_vcf
#SBATCH -t 148:00:00
#SBATCH -c 48
#SBATCH --mem=720G
#SBATCH --ntasks=1
#SBATCH --partition=bmh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

cd VCF_withKistler_noHerbarium

module load bcftools

bcftools annotate --rename-chrs ../Chr_renaming_map.txt combined_cassava_mar2025_final_withKistler_noHerbarium.vcf.gz -Oz -o combined_cassava_mar2025_final_withKistler_noHerbarium_chrrenamed.vcf.gz

echo "renaming finished"

module load jdk

java -Xmx8g -jar ~/snpEff/snpEff.jar -v Cassava_v8 combined_cassava_mar2025_final_withKistler_noHerbarium_chrrenamed.vcf.gz > Cassava_complete_ann_apr2025.vcf.gz
