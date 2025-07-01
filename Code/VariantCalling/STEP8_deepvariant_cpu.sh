#!/bin/bash -l
#SBATCH -J DeepVariant
#SBATCH -t 96:00:00
#SBATCH -c 48
#SBATCH --mem=360G
#SBATCH --ntasks=1
#SBATCH --partition=bmm
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kehzhao@ucdavis.edu
#SBATCH --output=./report/slurm-%j.out

READS=$1
echo "$READS"

module load apptainer

# Command-line arguments

ref_genome=./GCF_001659605.2_M.esculenta_v8_genomic.fasta

rm -r ./intermediate_results_dir/$READS
mkdir ./intermediate_results_dir/$READS

singularity exec /home/icanders/deepvariant_gpu.sif run_deepvariant \
   --model_type=WGS \
   --ref=$ref_genome \
   --reads=./sorted_bam/$READS.fix.markdup.bam \
   --output_vcf=./deepvariant_output/$READS.vcf.gz \
   --output_gvcf=./deepvariant_output/$READS.gvcf.gz \
   --intermediate_results_dir ./intermediate_results_dir/$READS \
   --num_shards=48 \

rm -r ./intermediate_results_dir/$READS
