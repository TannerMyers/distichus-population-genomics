#!/usr/bin/env bash

# This script runs the entire Stacks pipeline from ustacks through populations #

# Assign variables

## Set command line arguments for M and n -- the number of mismatches allowed 
## between stacks within individuals and mismatches allowed between stacks
## between individuals
M=$1
n=$2

SAMPLE_DIR=/scratch/phyletica/distichus/samples

POP_MAP=/scratch/phyletica/distichus/info/popmap.tsv

SAMPLES=`awk 'BEGIN {OFS = FS} {print $1}' $POP_MAP` 

OUT_DIR=/scratch/phyletica/distichus/stacks.denovo/uncleaned.M4m3N3

# Load 
module load stacks

id=1
for sample in $SAMPLES 
do
	ustacks -f $SAMPLE_DIR/${sample}.1.fq.gz -i $id --name $sample -o $OUT_DIR -M $M -m $n -N $n --disable-gapped -p 8
	let "id+=1"
done

cstacks -P $OUT_DIR \
	-M $POP_MAP \
	-n $n \
	-p 8	

sstacks -P $OUT_DIR \
	-M $POP_MAP \
	-p 8 

tsv2bam -P $OUT_DIR \
        -M $POP_MAP \
	--pe-reads-dir $SAMPLE_DIR \
        -t 8 

gstacks -P $OUT_DIR \
	-M $POP_MAP \
        -t 8 

#populations -P /scratch/phyletica/distichus/stacks.denovo/uncleaned \
#	-M /scratch/phyletica/distichus/info/popmap.tsv \
#	-t 8 \
#	-R 0.65 \
#	--min-mac 2 \
#	--write-single-snp \
#	--fasta-samples \
#	--fasta-loci \
#	--vcf \
#	--structure \
#	--hwe \
#	--fstats \
#	--phylip-var \
#	--plink
	 
