#!/usr/bin/env bash

# Load necessary modules
module load stacks
module load bwa/0.7.17 
module load samtools
module load python/3.6.4

BWA_DB=/scratch/phyletica/distichus/genome/bwa/

CAT_LOCI=catalog.fa.gz 

STACKS_DIR=/scratch/phyletica/distichus/stacks.denovo/uncleaned.M4n3

OUT_DIR=$STACKS_DIR/integrated-alignments/

    bwa mem

    stacks-integrate-alignments -P STACKS_DIR \
        -B BAM_FILE \
        -O OUT_DIR    