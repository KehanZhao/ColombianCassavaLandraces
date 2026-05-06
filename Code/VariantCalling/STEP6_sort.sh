#!/bin/bash -l
#SBATCH -J sort
#SBATCH -t 96:00:00
#SBATCH --mem=500G
#SBATCH --ntasks=32
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

module load samtools

cd bam

x=$1
echo "$x"

samtools sort -n -@5 -m 200G -o ./sorted_bam/$x.sorted.bam $x.bam 
