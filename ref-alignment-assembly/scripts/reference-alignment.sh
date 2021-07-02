#!/usr/bin/env bash

################################################################

# This script uses BWA-MEM to align the demultiplexed ddRADseq 
# fastq files from process_radtags to the Anolis carolinensis
# reference genome indexed by the script "bwa_index.sh". This 
# will provide an alternative to Stacks for RAD loci assembly.

################################################################

# Designate sample ID variables

    # Population map file containing individual IDs that for loop will iterate over.
    POP_MAP=/scratch/phyletica/distichus/info/popmap.tsv


    # Isolate individual IDs in population map to use for looping variable
    SAMPLES=`awk 'BEGIN {OFS = FS} {print $1}' $POP_MAP`

# Name input & output variables

    # Directory holding reference genome database for Anolis carolinensis
    BWA_DB=/scratch/phyletica/distichus/genome/bwa/anocar

    # Directory containing demultiplexed fastq files from process_radtags
    SAMPLE_DIR=/scratch/phyletica/distichus/samples/

    # Directory to output files produced by bwa mem
    OUT_DIR=/scratch/phyletica/distichus/alignments 


for sample in $SAMPLES; 
    do
        qsub -N $sample.refalign.test -d /scratch/phyletica/distichus/scripts -q gen28 -W group_list=jro0014_lab -W x=FLAGS:ADVRES:jro0014_s28 -l nodes=1:ppn=8,mem=80gb,walltime=12:00:00 <<<" 

            # Load conda environment with bwa
            source ~/mambaforge/etc/profile.d/conda.sh # Or path to where your conda is
            conda activate genomics_env
    	    	

            # Provide forward and reverse reads
            READ1=$SAMPLE_DIR$sample.1.fq.gz
            READ2=$SAMPLE_DIR$sample.2.fq.gz

	    # Align paired read fastq files to the indexed Anolis carolinensis reference, obtaining individual .sam files 
            #bwa mem -t 8 /scratch/phyletica/distichus/genome/bwa/anocar $SAMPLE_DIR$sample.1.fq.gz $SAMPLE_DIR$sample.2.fq.gz -o $OUT_DIR/bwa-outputs/$sample.sam
	    ## Should be bwa mem -t 8 $BWA_DB $READ1 $READ2 -o $OUT_DIR/bwa-outputs/$sample.sam

	    # Use samtools to convert to .bam, sort, and index
	    samtools view -b $OUT_DIR/bwa-outputs/$sample.sam -o $OUT_DIR/bam-files/$sample.bam 
	    samtools sort -o $OUT_DIR/bam-files/$sample.sorted.bam $OUT_DIR/bam-files/$sample.bam
	    samtools index $OUT_DIR/bam-files/$sample.sorted.bam

	    "
    done
