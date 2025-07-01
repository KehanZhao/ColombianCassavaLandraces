#!/bin/bash -l
#SBATCH -J merge
#SBATCH -t 96:00:00
#SBATCH --mem=500G
#SBATCH --ntasks=32
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

# This is examples of merging multiple bams of the same biosample
module load samtools

cd bam

samtools merge  -o SRR2847405_merged.bam SRR2847406.bam SRR2847405.bam
samtools merge  -o SRR2847384_merged.bam SRR2847384.bam SRR2847387.bam
samtools merge  -o SRR486604_merged.bam SRR486604.bam SRR486605.bam
samtools merge  -o SRR1261871_merged.bam SRR2847442.bam SRR2847441.bam SRR1261871.bam
samtools merge  -o SRR1261886_merged.bam SRR2847450.bam SRR2847449.bam SRR1261886.bam
samtools merge  -o SRR1261873_merged.bam SRR2847448.bam SRR2847447.bam SRR1261873.bam
samtools merge  -o SRR1261872_merged.bam SRR1261872.bam SRR2847443.bam SRR2847444.bam SRR2847445.bam SRR2847446.bam
samtools merge  -o SRR1261917_merged.bam SRR2847452.bam SRR2847451.bam SRR1261917.bam
samtools merge  -o SRR1261867_merged.bam SRR2847415.bam SRR1261867.bam
samtools merge  -o SRR2847393_merged.bam SRR2847394.bam SRR2847393.bam
samtools merge  -o SRR2847380_merged.bam SRR2847381.bam SRR2847380.bam
samtools merge  -o SRR1261958_merged.bam SRR2847379.bam SRR1261958.bam
samtools merge  -o SRR2847474_merged.bam SRR2847474.bam SRR2847473.bam
samtools merge  -o SRR2847468_merged.bam SRR2847468.bam SRR2847467.bam
samtools merge  -o SRR2847459_merged.bam SRR2847459.bam SRR2847458.bam
samtools merge  -o SRR2847396_merged.bam SRR2847396.bam SRR2847395.bam
samtools merge  -o SRR2847454_merged.bam SRR2847454.bam SRR2847453.bam
samtools merge  -o SRR2847456_merged.bam SRR2847456.bam SRR2847455.bam
samtools merge  -o SRR2847440_merged.bam SRR2847440.bam SRR2847439.bam
samtools merge  -o SRR2847429_merged.bam SRR2847429.bam SRR2847428.bam SRR2847427.bam
samtools merge  -o SRR2847432_merged.bam SRR2847432.bam SRR2847431.bam SRR2847430.bam
samtools merge  -o SRR2847435_merged.bam SRR2847435.bam SRR2847434.bam SRR2847433.bam
samtools merge  -o SRR2847438_merged.bam SRR2847438.bam SRR2847437.bam SRR2847436.bam
samtools merge  -o SRR2847466_merged.bam SRR2847466.bam SRR2847465.bam
samtools merge  -o SRR2847392_merged.bam SRR2847392.bam SRR1261865.bam
samtools merge  -o SRR2847472_merged.bam SRR2847472.bam SRR2847471.bam SRR2847470.bam SRR2847469.bam SRR1261959.bam
samtools merge  -o E150016662_L01_104_merged.bam E150016662_L01_104.bam E150016689_L01_104.bam
samtools merge  -o E150016662_L01_122_merged.bam E150016662_L01_122.bam E150016689_L01_122.bam
samtools merge  -o E150016662_L01_123_merged.bam E150016662_L01_123.bam E150016689_L01_123.bam
samtools merge  -o E150016662_L01_102_merged.bam E150016662_L01_102.bam E150016689_L01_102.bam

