#!/usr/bin/env bash

# This script runs the entire Stacks pipeline from ustacks through populations #

# Assign variables

## Set command line arguments for M and n -- the number of mismatches allowed 
## between stacks within individuals and mismatches allowed between stacks
## between individuals
#M=$1
#n=$2 
#POP_MAP=$3
#OUT_DIR=$4
#POPULATION_DIR=$5
POP_MAP=$1
OUT_DIR=$2
POPULATIONS_DIR=$3

#SAMPLE_DIR=/scratch/phyletica/distichus/samples

#SAMPLES=`awk 'BEGIN {OFS = FS} {print $1}' $POP_MAP`

#BWA_DB=/scratch/phyletica/distichus/genome/bwa/anocar

# Load necessary modules 
module load stacks
module load python
module load samtools
module load bwa

#id=1
#for sample in $SAMPLES 
#do
#	ustacks -f $SAMPLE_DIR/${sample}.1.fq.gz -i $id --name $sample -o $OUT_DIR -M $M -m $n -N $n --disable-gapped -p 8
#	let "id+=1"
#done

#cstacks -P $OUT_DIR \
#	-M $POP_MAP \
#	-n $n \
#	-p 8	

#sstacks -P $OUT_DIR \
#	-M $POP_MAP \
#	-p 8 

#tsv2bam -P $OUT_DIR \
#        -M $POP_MAP \
#		--pe-reads-dir $SAMPLE_DIR \
#        -t 8 

#gstacks -P $OUT_DIR \
#	-M $POP_MAP \
#        -t 8 

#bwa mem -t 8 $BWA_DB $OUT_DIR/catalog.fa.gz |
#    samtools view -b |
#    samtools sort --threads 4 > $OUT_DIR/aligned_catalog.bam

# Confused about what "stacks_dir" refers to with the "-P" flag
#stacks-integrate-alignments -P $OUT_DIR \
#    -B $OUT_DIR/aligned_catalog.bam \
#    -O $OUT_DIR/integrated-alignment

populations -P $OUT_DIR \
	-O $POPULATIONS_DIR \
	-M $POP_MAP \
	-t 10 \
	-r 0.3 \
	--min-mac 2 --write-single-snp \
	--phylip-var \
	--plink \
	--structure
#       -R 0.65 \
#       --write-single-snp \
#       --fasta-samples \
#       --fasta-loci \
#	--hwe \
#	--fstats \
#	--vcf	 
