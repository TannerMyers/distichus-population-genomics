#!/usr/bin/env bash

module load stacks


POP_STACKS_DIR=/scratch/phyletica/distichus/stacks.denovo/population-stacks.denovo/

POP_MAP_DIR=/scratch/phyletica/distichus/info/lineage-popmaps

POPULATIONS="brevirostris
	distichus
	dominicensis
	favillarum
	ignigularis
	ravitergum
	Tiburon"

for pop in $POPULATIONS
	do 
		cd $POP_STACKS_DIR$pop
		populations -P $POP_STACKS_DIR$pop -M $POP_MAP_DIR/popmap_${pop}.tsv -t 8 --vcf
	done
