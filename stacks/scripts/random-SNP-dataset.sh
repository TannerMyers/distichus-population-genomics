#!/usr/bin/env bash

# Load modules
module load stacks

# Name variables
OUT_DIR=$1
POPULATIONS_DIR=$2
POP_MAP=$3



[[ -d $POPULATIONS_DIR/run_%j ]] || mkdir -p $POPULATIONS_DIR/run_%j 
# the %j is slurm for the job number -- figure out how to do this with PBS
# I'm just going to use 1-100 and use that for iterating my for loop and naming temporary directories to contain outputs


for ; do

qsub

	populations --in-path $OUT_DIR \
		--out-path $POPULATIONS_DIR \
		--popmap $POP_MAP \
		--threads \
		--ordered-export \
		--write-random-snp \
		--min-mac 2 \
		--vcf --plink 

	cp $POPULATIONS_DIR/populations.structure $structure_folder/populations.structure.%j

done
