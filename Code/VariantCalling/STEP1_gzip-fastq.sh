#!/bin/bash -l
#SBATCH -J gzip
#SBATCH -t 96:00:00
#SBATCH --mem=8G 
#SBATCH --ntasks=1
#SBATCH --partition=bmm
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

SRR=$1 # The accession number of the reads to be gzipped
cd fastq # Change to directory where fastq files are located
# First, check if the files are already partially gzipped (avoids requeuing errors)
# Delete the partially gzipped file if it exists
if [ -e ${SRR}_1.fastq.gz ]
then
  rm ${SRR}_1.fastq.gz
fi
if [ -e ${SRR}_2.fastq.gz ]
then
  rm ${SRR}_2.fastq.gz
fi
cd .. # Change back to directory containing this file (gzip-fastq.sh)
# gzip the forward and the reverse reads
gzip ./fastq/${SRR}_1.fastq
gzip ./fastq/${SRR}_2.fastq
