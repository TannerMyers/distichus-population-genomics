#! /usr/bin/env bash

#################################################################
# This script assesses the mapping quality of .bam files output #
# by "reference-alignment.sh" with the Java tool Qualimap       #
#################################################################
# Load conda environment with bwa, samtools, and qualimap
    source ~/mambaforge/etc/profile.d/conda.sh # Or path to where your conda is
    conda activate genomics_env
    	    	
qualimap multi-bamqc --run-bamqc --data /scratch/phyletica/distichus/alignments/qualimap-input-information.tsv -gff /scratch/phyletica/distichus/genome/gtf/Anolis_carolinensis.AnoCar2.0v2.104.gtqualimap multi-bamqc --run-bamqc --data /scratch/phyletica/distichus/alignments/qualimap-input-information.tsv -gff /scratch/phyletica/distichus/genome/gtf/Anolis_carolinensis.AnoCar2.0v2.104.gtff