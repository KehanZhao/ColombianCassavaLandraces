#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=1008g
#SBATCH --time=5-20:00:00
#SBATCH --job-name=beagle
#SBATCH --output=jobname.out.%j
#SBATCH --partition=bmh
eval "$(conda shell.bash hook)"
conda activate BEAGLE
Genome=$1

beagle -Xmx1000g gt=Cassava_complete_ann_7_2023_All.vcf out=Cassava_complete_ann_7_2023_All_beagle5 ne=1000  nthreads=100 window=50
#beagle -Xmx1000G gt=test5.vcf out=Cassava_complete_ann_7_2023_beagle5
