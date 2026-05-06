#!/bin/bash -l
#SBATCH -J fix
#SBATCH -t 96:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=5
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

SRR=$1 # The accession number of the mapped reads to be fixed
echo "$SRR"

cd sorted_bam
module load samtools

# Fix reads
samtools fixmate -m $SRR.sorted.bam - | samtools sort -@5 -o $SRR.fix.bam
# Mark duplicates
samtools markdup -@5 -s $SRR.fix.bam - | samtools sort -@5 -o $SRR.fix.markdup.bam
# Generate index file
samtools index $SRR.fix.markdup.bam
