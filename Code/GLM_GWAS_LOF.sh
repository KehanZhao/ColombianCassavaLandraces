#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=88g
#SBATCH --time=8:00:00
#SBATCH --job-name=Tassel
#SBATCH --output=jobname.out.%j
#SBATCH --partition=bmh
eval "$(conda shell.bash hook)"
conda activate Landrace_tools

run_pipeline.pl  -Xms80g  -Xmx80g -fork1 -importGuess $1  \
-fork2 -importGuess $2 \
-fork3 -importGuess Colombia_Batch1234_PC.txt \
-fork4 -importGuess Colombia_Batch1234_kinship.txt \
-combine5 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin -endplugin -export LOF_glm_output
##Kinship is not utilized in this model, even though the file is loaded

