#!/bin/bash -l
#SBATCH -J fasterq-dump
#SBATCH -t 96:00:00
#SBATCH --mem=20G 
#SBATCH --ntasks=1
#SBATCH --partition=bmm
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

SRR=$1 # The accession number to be downloaded
echo "$SRR"

# Checks if the fastq files for this SRR already exist and deletes them if so.
# Avoids requeing errors where the incomplete file is not deleted before
# the program is run, which can happen on bmm
if [ -e ./fastq/${SRR}_1.fastq -o -e ./fastq/${SRR}_2.fastq -o -e ${SRR}.fastq ]
then
  rm ./fastq/${SRR}*.fastq
fi

module load sratoolkit
fasterq-dump --split-files $SRR --outdir ./fastq

# Once both forward and reverse reads are downloaded,
# the files can be forwarded to another sbatch file to be gzipped.
# Make sure gzip-fastq exists before running, or delete
sbatch STEP1_gzip-fastq.sh $SRR
