#!/usr/bin/env bash

# Load conda environment containing Admixture v. 1.3.0
source ~/miniconda3/etc/profile.d/conda.sh # Or path to where your conda is
conda activate population-genomics

# Name variables
INPUT_FILE=$1

# Run Admixture on a number of "K" value.
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15; 
	do
	    # Run admixture with cross validation
	    admixture --cv $INPUT_FILE $K | tee log${K}.out

	done
