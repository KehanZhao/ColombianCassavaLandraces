#!/bin/bash -l
#SBATCH -J admixture_escu_noHu
#SBATCH -t 720:00:00
#SBATCH --mem=200G
#SBATCH --ntasks=2
#SBATCH --partition=bmh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

PROJECT=/home/kehan/project/Cassava/Ramu_extra
VCF_DIR=$PROJECT/VCF_withKistler_noHerbarium

# filter snps
vcftools --gzvcf $VCF_DIR/Cassava_final_1152accessions_apr2025.vcf.gz \
  --minQ 30 \
  --positions $PROJECT/VCF/XGboost_Techrep_SitePASS.txt \
  --maf 0.05 --stdout --recode --recode-INFO-all \
  > $VCF_DIR/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost.vcf

echo "snp filter finished"

# filter accessions
conda activate bcftools

bcftools view \
  -S $VCF_DIR/Accession_list_filtering/Esculenta_noHu.txt \
  $VCF_DIR/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost.vcf \
  -Ov -o $VCF_DIR/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu.vcf

conda deactivate

echo "accession filter finished"

# make bed
conda activate Plink

plink \
  --vcf $VCF_DIR/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu.vcf \
  --chr 1-18 --make-bed --double-id --vcf-half-call missing --allow-extra-chr \
  --out $VCF_DIR/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu

plink \
  --bfile $VCF_DIR/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu \
  --geno 0.1 --make-bed \
  --out $VCF_DIR/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct

conda deactivate

echo "make bed completed"

conda activate ADMIXTURE

for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30; \
do admixture --cv $VCF_DIR/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct.bed $K | tee Admixture_output/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost_esculenta_excludingHu_geno10pct_log${K}.out; done
