#!/bin/bash -l
#SBATCH -J trimmomatic
#SBATCH -t 96:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

eval "$(conda shell.bash hook)"
conda activate trimmomatic

cd fastq/ # Change to directory where fastq files are located

SRR=$1 # The accession number of the reads to be trimmed
echo "$SRR"

# Specify phred33 or phred64 based on encoding if known. This prevents reader error and improves speed
trimmomatic PE -threads 4 -phred33 \
  ${SRR}_1.fastq.gz ${SRR}_2.fastq.gz \
  ${SRR}_1.trimmed.fastq.gz ${SRR}_1un.trimmed.fastq.gz \
  ${SRR}_2.trimmed.fastq.gz ${SRR}_2un.trimmed.fastq.gz \
  SLIDINGWINDOW:4:20

# Use the trimmed PAIRED reads for next steps
