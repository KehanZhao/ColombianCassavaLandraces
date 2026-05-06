#!/bin/bash -l
#SBATCH -J Plink_LD
#SBATCH -t 124:00:00
#SBATCH --mem=200G
#SBATCH --ntasks=1
#SBATCH --partition=bmh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

your_input=$1

VCF_DIR=/home/kehan/project/Cassava/Ramu_extra/VCF_withKistler_noHerbarium

# filter accessions
conda activate bcftools

bcftools view \
  -S $VCF_DIR/Accession_list_filtering/${your_input}.txt \
  $VCF_DIR/Cassava_final_1152accessions_apr2025_Qual30_maf5pct_xgboost.vcf \
  -Ov -o $VCF_DIR/${your_input}_accesions_filtered.vcf

conda deactivate

conda activate Plink

plink \
  --vcf ../../VCF_withKistler_noHerbarium/${your_input}_accesions_filtered.vcf \
  --make-bed --double-id --vcf-half-call missing \
  --out ../../VCF_withKistler_noHerbarium/${your_input}_accesions_filtered

plink \
  --bfile ../../VCF_withKistler_noHerbarium/${your_input}_accesions_filtered \
  --r2 --ld-window-r2 0 --ld-window 99999 --ld-window-kb 100 --thin-count 50000 \
  --out ${your_input}_accesions_filtered_window100kb_thincount50000_may2025
