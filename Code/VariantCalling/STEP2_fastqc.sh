#!/bin/bash -l
#SBATCH -J fastqc
#SBATCH -t 96:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --partition=bmm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

SRR=$1 # The accession number of the reads to generate reports of

# First, check if fatqc files already exist to avoid requeing errors

# Make a list to easily loop through and check if they exist
fastqc_files=(${SRR}_1_fastqc.html ${SRR}_1_fastqc.zip ${SRR}_2_fastqc.html ${SRR}_2_fastqc.zip)

# Change to directory where fastq files are located
cd fastq

# Loop through list and delete the fastqc files that already exist
for FILE in "${fastqc_files[*]}"
do
  if [ -e $FILE ]
  then
    rm $FILE
  fi
done

# Generate fatqc files for forward and reverse reads
module load fastqc/0.11.9
zcat ${SRR}_1.fastq.gz | fastqc stdin:${SRR}_1
zcat ${SRR}_2.fastq.gz | fastqc stdin:${SRR}_2
