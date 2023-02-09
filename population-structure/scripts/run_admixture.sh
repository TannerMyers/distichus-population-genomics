#!/usr/bin/env bash

# Load conda environment containing Admixture v. 1.3.0
source ~/mambaforge/etc/profile.d/conda.sh # Or path to where your conda is
conda activate popgen_env

# Name variables
INPUT_FILE=$1

CV_LOG=$2

# Run Admixture on a number of "K" value.
for K in 2 3 4 5 6 7 8 9 10; 
	do
	    # Run admixture with cross validation
	    admixture --cv $INPUT_FILE $K | tee log${K}.out >> $CV_LOG

	done
