#!/bin/bash -l
#SBATCH -J vcf_filter_esculenta
#SBATCH -t 48:00:00
#SBATCH --mem=500G
#SBATCH --ntasks=1
#SBATCH --partition=bmh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=../../report/slurm-%j.out

conda activate bcftools

bcftools view -S ../../VCF_withKistler_noHerbarium/esculenta_lines_noHerbarium.txt ../../VCF_withKistler_noHerbarium/All_accessions_noHerbarium_filtered_Qual30_maf5pct_xgboost.vcf -Ov -o ../../VCF_withKistler_noHerbarium/Esculenta_accesions_filtered.vcf

conda deactivate

echo "Filter finished"

conda activate Plink
plink --vcf /home/kehan/project/Cassava/Ramu_extra/VCF_withKistler_noHerbarium/Esculenta_accesions_filtered.vcf --chr 1-18 --make-bed --double-id --vcf-half-call missing --allow-extra-chr --out /home/kehan/project/Cassava/Ramu_extra/VCF_withKistler_noHerbarium/Esculenta_accesions_filtered

# PCA
plink --bfile /home/kehan/project/Cassava/Ramu_extra/VCF_withKistler_noHerbarium/Esculenta_accesions_filtered --pca --out ./PCA_output/Esculenta_accesions_withKistler_noHerbarium_filtered