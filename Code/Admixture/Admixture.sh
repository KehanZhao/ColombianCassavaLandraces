#!/bin/bash -l
#SBATCH -J admixture_escu_noHu
#SBATCH -t 720:00:00
#SBATCH --mem=200G
#SBATCH --ntasks=2
#SBATCH --partition=bmh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

#conda activate bcftools

#bcftools view -S ../../VCF/Esculenta_without_Hu_list.txt ../../VCF/All_accessions_except_Hu_filtered_GQ20_Qual30_maf5pct_xgboost.vcf -Ov -o ../../VCF/Esculenta_accesions_without_Hu_filtered.vcf

#conda deactivate

#echo "Filter finished"

# make bed
conda activate Plink

#plink --vcf /home/kehan/project/Cassava/Ramu_extra/VCF/Esculenta_accesions_without_Hu_filtered.vcf --chr 1-18 --make-bed --double-id --vcf-half-call missing --allow-extra-chr --out /home/kehan/project/Cassava/Ramu_extra/VCF/Esculenta_accesions_without_Hu_filtered

plink --bfile /home/kehan/project/Cassava/Ramu_extra/VCF/Esculenta_accesions_without_Hu_filtered --geno 0.1 --make-bed --out /home/kehan/project/Cassava/Ramu_extra/VCF/Esculenta_accesions_without_Hu_filtered_geno10pct

conda deactivate

echo "make bed completed"

conda activate ADMIXTURE

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30; \
do admixture --cv /home/kehan/project/Cassava/Ramu_extra/VCF/Esculenta_accesions_without_Hu_filtered_geno10pct.bed $K | tee Admixture_output/Esculenta_accesions_without_Hu_filtered_log${K}.out; done
