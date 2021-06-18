#!/usr/bin/env bash

################################################################

# This script uses BWA-MEM to align the demultiplexed ddRADseq 
# fastq files from process_radtags to the Anolis carolinensis
# reference genome indexed by the script "bwa_index.sh". This 
# will provide an alternative to Stacks for RAD loci assembly.

################################################################

# Name variables
    # Directory holding reference genome database for Anolis carolinensis
    BWA_DB=/scratch/phyletica/distichus/genome/bwa/ 

    # Directory containing demultiplexed fastq files from process_radtags
    SAMPLE_DIR=/scratch/phyletica/distichus/samples/

    # Directory to output files produced by bwa mem
    OUT_DIR=/scratch/phyletica/distichus/alignments

    # Population map file containing individual IDs that for loop will iterate over.
    POP_MAP=/scratch/phyletica/distichus/info/popmap.tsv

    # Isolate individual IDs in population map to use for looping variable
    SAMPLES=`awk 'BEGIN {OFS = FS} {print $1}' $POP_MAP`


for sample in $SAMPLES; 
    do
        qsub -N $sample.bwa \
            -d /scratch/phyletica/distichus/scripts \
            -q gen28 \
		    -W group_list=jro0014_lab \
            -W x=FLAGS:ADVRES:jro0014_s28 \
    		-l nodes=1:ppn=8,mem=80gb,walltime=20:00:00 <<<"            

            # Load conda environment with bwa
            source ~/miniconda3/etc/profile.d/conda.sh # Or path to where your conda is
            conda activate genomics_env
    
            # Provide forward and reverse reads
            READ1=$SAMPLE_DIR$sample.1.fq.gz
            READ2=$SAMPLE_DIR$sample.2.fq.gz

            bwa mem -t 8 $BWA_DB \
                $READ1 $READ2 \
                2> $OUT_DIR/logs/$sample-bwa.err \
                > $OUT_DIR/bwa-outputs/$sample.sam
            "
    done