#!/usr/bin/env bash

# Name variables
INPUT_DIR=/scratch/phyletica/distichus/stacks.denovo/uncleaned.M4m3N3/
POP_MAP=/scratch/phyletica/distichus/info/popmap_cleaned_MacGuigan_no-brev_inds-as-pops.tsv
POPULATIONS_DIR=/scratch/phyletica/distichus/population-structure-datasets

for ((i=1; i<=100;i++)); 
	do

		qsub -N admixture-dataset$i \
			-n -d /scratch/phyletica/distichus/scripts \
			-q gen28 \
			-W group_list=jro0014_lab \
			-W x=FLAGS:ADVRES:jro0014_s28 \
			-l nodes=1:ppn=10,mem=100gb,walltime=24:00:00 <<<"

				# Load modules
				module load stacks

				# Make a directory for this populations run's number if it doesn't already exist 
				[[ -d $POPULATIONS_DIR/Run$i ]] || mkdir -p $POPULATIONS_DIR/Run$i 

				# Run populations to get plink (.ped), structure, and .vcf files
				populations --in-path $INPUT_DIR \
					--out-path $POPULATIONS_DIR/Run$i \
					--popmap $POP_MAP \
					--threads 10 \
					-R 0.7 \
					--ordered-export \
					--write-random-snp \
					--min-mac 2 \
					--vcf --plink --structure

				# Copy populations output files to directory that will contain the populations
				# output files for all datasets
				cp $POPULATIONS_DIR/Run$i/populations.snps.vcf $POPULATIONS_DIR/run$i.populations.snps.vcf
				cp $POPULATIONS_DIR/Run$i/populations.plink.ped $POPULATIONS_DIR/run$i.populations.plink.ped
				cp $POPULATIONS_DIR/Run$i/populations.plink.map $POPULATIONS_DIR/run$i.populations.plink.map
				cp $POPULATIONS_DIR/Run$i/populations.structure $POPULATIONS_DIR/run$i.populations.structure
				cp $POPULATIONS_DIR/Run$i/populations.log $POPULATIONS_DIR/run$i.populations.log
				"
	done
