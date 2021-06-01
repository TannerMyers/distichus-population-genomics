#!/usr/bin/env bash

# Submit R script as job

	qsub -N conStruct-trial3 \
		-n \
		-d /home/tcm0036/distichus/population-structure/conStruct \
		-q gen28 \
		-W group_list=jro0014_lab \
		-W x=FLAGS:ADVRES:jro0014_s28 \
		-l nodes=1:ppn=10,mem=80gb,walltime=72:00:00 <<<"
	
	# Load conda environment with R
	source ~/miniconda3/etc/profile.d/conda.sh # Or path to where your conda is
	conda activate R_env

	~/miniconda3/envs/R_env/bin/Rscript conStruct_analysis.R
	"
