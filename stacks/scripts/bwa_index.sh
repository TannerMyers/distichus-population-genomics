#!/usr/bin/env bash

# Download genome from Ensembl FTP site beforehand
genome_fa=/scratch/phyletica/distichus/genome/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa.gz

# Load bwa module
module load bwa/0.7.17  

# Create index for reference genome
bwa index -p /scratch/phyletica/distichus/genome/bwa/anocar $genome_fa
