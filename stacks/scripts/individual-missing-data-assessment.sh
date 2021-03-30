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

# Load modules
module load vcftools 

# Loop through directories containing Stacks outputs for individual 
# populations to run vcftools 
for pop in ; # Need to figure out how to iterate through directories list
    do 
        cd $POP_DIR$pop
        vcftools \ 
            --vcf populations.snps.vcf \ 
            --missing-indv 
    done
