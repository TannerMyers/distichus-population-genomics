#!/usr/bin/env bash

################################################################

# This script uses BWA-MEM to align the demultiplexed ddRADseq 
# fastq files from process_radtags to the Anolis carolinensis
# reference genome indexed by the script "bwa_index.sh". This 
# will provide an alternative to Stacks for RAD loci assembly.

#Note: Prior to running this script, users should create a 
# conda environment with the packages involved here if they 
# are not using their HPC's modules

################################################################

set -e

# Population map file containing individual IDs that for loop will iterate over.
popmap=/home/tcm0036/distichus-ddRAD/info/distichus-popmap.tsv

	# Isolate individual IDs in population map to use as looping variables
	samples=`awk 'BEGIN {OFS = FS} {print $1}' $popmap`

for ID in $samples;
	do
	sbatch --job-name=$ID-alignment --partition=jro0014_amd --cpus-per-task=8 --mem=80G --time=20-24:00:00 --wrap="	

	# Load modules 
	module load bwa
	module load samtools

	# Name variables 
	# Path to Anolis carolinensis genome fasta
	genome=/home/tcm0036/distichus-ddRAD/genome/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa
	# Path to including prefix of database generated by bwa index
	bwa_db=/home/tcm0036/distichus-ddRAD/genome/bwa/anocar

	fq1=/home/tcm0036/distichus-ddRAD/samples/$ID.1.fq.gz
	fq2=/home/tcm0036/distichus-ddRAD/samples/$ID.2.fq.gz
	sam=/scratch/tcm0036/distichus-ddRAD/alignments/results/sam/$ID.aligned.sam
	bam=/scratch/tcm0036/distichus-ddRAD/alignments/results/bam/$ID.aligned.bam
	sorted_bam=/scratch/tcm0036/distichus-ddRAD/alignments/results/bam/$ID.aligned.sorted.bam

	bwa mem -t 8 \$bwa_db \$fq1 \$fq2 > \$sam
	samtools view -b \$sam > \$bam
	samtools sort -o \$sorted_bam \$bam
	samtools index \$sorted_bam
	"

done
