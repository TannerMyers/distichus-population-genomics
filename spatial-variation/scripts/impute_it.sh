#!/usr/bin/env bash
#SBATCH --job-name linkImputeR
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tcm0036@auburn.edu
#SBATCH --time=48:00:00
#SBATCH --mem 75G
#SBATCH --partition bac0071_amd



# Load conda environment with bcftools, vcftools, and vcflib
source ~/mambaforge/etc/profile.d/conda.sh # Or path to where your conda is
conda activate genomics_env

java -jar LinkImputeR/LinkImputeR.jar -s accuracy.ini
