#!usr/bin/env bash

# Load conda environment with R
source ~/miniconda3/etc/profile.d/conda.sh # Or path to where your conda is
conda activate R_env

# Submit R script as job
Rscript conStruct_analysis.R