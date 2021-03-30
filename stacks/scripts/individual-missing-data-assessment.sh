#! /usr/bin/env bash

##################################################################
# This script runs vcftools with the `--missing-indv` flag on    #
# the `populations.snps.vcf` to calculate the percentage of data #
# that is missing at an individual level. Then, I calculate the  #
# mean % missing data and list the individuals with missing data #
# higher than the mean for the population level.                 #
# This approach follows the protocol of Cerca et al. 2021.       #
##################################################################

# Name variables 
POP_DIR=/scratch/phyletica/distichus/stacks.denovo/population-stacks.denovo

POPULATIONS=`ls $POP_DIR`

# Load modules
module load vcftools/v0.1.17          

# Loop through directories containing Stacks outputs for individual 
# populations to run vcftools 
for pop in $POPULATIONS; 
	do 
		# Move into directory containing Stacks output files for a given population.	
		cd $POP_DIR/$pop
	
		# Estimate % missing data for each individual in the population.
		# This results in two output files: `out.imiss` and `out.log`
        #vcftools --vcf populations.snps.vcf --missing-indv 

		# Extract individual IDs and % missing data to single file
		#awk '{print $1, $5}' out.imiss >> $POP_DIR/missing_data_total

	done
