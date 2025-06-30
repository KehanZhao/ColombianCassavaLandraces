#!/bin/bash -l
#SBATCH -J vcf_filter
#SBATCH -t 48:00:00
#SBATCH --mem=500G
#SBATCH --ntasks=1
#SBATCH --partition=bmh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=../../report/slurm-%j.out


vcftools --gzvcf /home/kehan/project/Cassava/Ramu_extra/VCF_withKistler_noHerbarium/combined_cassava_mar2025_final_withKistler_noHerbarium_chrrenamed.vcf.gz --minQ 30 --positions /home/kehan/project/Cassava/Ramu_extra/VCF/XGboost_Techrep_SitePASS.txt --maf 0.05 --stdout --recode --recode-INFO-all  > /home/kehan/project/Cassava/Ramu_extra/VCF_withKistler_noHerbarium/All_accessions_noHerbarium_filtered_Qual30_maf5pct_xgboost.vcf

echo "Filter finished"

conda activate Plink
plink --vcf /home/kehan/project/Cassava/Ramu_extra/VCF_withKistler_noHerbarium/All_accessions_noHerbarium_filtered_Qual30_maf5pct_xgboost.vcf --chr 1-18 --make-bed --double-id --vcf-half-call missing --allow-extra-chr --out /home/kehan/project/Cassava/Ramu_extra/VCF_withKistler_noHerbarium/All_accessions_noHerbarium_filtered_Qual30_maf5pct_xgboost

# PCA
plink --bfile /home/kehan/project/Cassava/Ramu_extra/VCF_withKistler_noHerbarium/All_accessions_noHerbarium_filtered_Qual30_maf5pct_xgboost --pca --out ./PCA_output/All_accessions_withKistler_noHerbarium_filtered_Qual30_maf5pct_xgboost