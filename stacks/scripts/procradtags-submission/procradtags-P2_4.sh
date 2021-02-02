#!/usr/bin/env bash

# Name variables
 ## Read variables are command line arguments. Enter your file names when you run script
READ1=$1
READ2=$2 
BARCODES=../info/barcodes_lane2_4.tsv
ENZYME1=sbfI 
ENZYME2=mspI 
OUT_DIR=/scratch/phyletica/distichus/samples/P2/
JOB_NAME=P2_4_processradtags

	qsub \
		-N $JOB_NAME \
                -q general \
                -W group_list=jro0014_lab \
                -W x=FLAGS:ADVRES:jro0014_lab \
                -d /scratch/phyletica/distichus/samples/ \
                 -l walltime=48:00:00 << HereDoc

		module load stacks
	
		process_radtags \
			-1 $READ1 \
			-2 $READ2 \
			-b $BARCODES \
			-o $OUT_DIR \
			--renz_1 $ENZYME1 \
			--renz_2 $ENZYME2 \
			-i gzfastq \
			--rescue \
			--quality \
			--clean \
			--inline_null &> process_radtags.lane2_4.oe 

HereDoc
