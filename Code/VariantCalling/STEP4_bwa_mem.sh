#!/bin/bash -l
#SBATCH -J bwa-mem
#SBATCH -t 96:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=32
#SBATCH --partition=bmm
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

eval "$(conda shell.bash hook)"
conda activate map_reads

cd fastq/

SRR=$1 # The accession number of the reads to be aligned
echo "$SRR"

# A variable can hold the path to the reference genome fasta file
REF=./GCF_001659605.2_M.esculenta_v8_genomic.fasta

# Call bwa-mem, and pipe output to samtools. Make output to the bam directory
bwa mem -M -t 32 -R "@RG\tID:$SRR\tSM:$SRR\tPL:ILLUMINA" $REF ${SRR}_1.trimmed.fastq.gz ${SRR}_2.trimmed.fastq.gz | samtools sort -n -@5 -o ../bam/$SRR.bam
