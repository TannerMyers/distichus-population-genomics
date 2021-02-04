#!/usr/bin/env bash 

M=$1

# Name --popmap variable
POP_MAP=/scratch/phyletica/distichus/info/popmap.test_samples.tsv

# Name --samples variable
SAMPLES=/scratch/phyletica/distichus/samples 

# Provide job submission name
JOB_NAME=denovomaptest

# Make an array for the values of M to test
#M=( 1 2 3 4 5 6 7 8 9 )

	# iterate over the array, M
#	for NUMBER in ${M[@]}; 
#			do
				# Name variable for output files 
		       		OUT_DIR=/scratch/phyletica/distichus/tests.denovo/stacks.M$M 
				# Move into directory for output files
				cd $OUT_DIR

			qsub -N $JOB_NAME$M \
                         -q general \
                         -W group_list=jro0014_lab \
                         -W x=FLAGS:ADVRES:jro0014_lab \
                         -d /scratch/phyletica/distichus/samples/ \
                         -l walltime=48:00:00 << HereDoc
				
				module load stacks

				# run denovo_map for those samples for the value of M being tested
						denovo_map.pl \
							-M $M \
							-n $M \
							--out-path $OUT_DIR \
							--popmap $POP_MAP \
							--samples $SAMPLES \
							--paired &> denovomap.oe 	

HereDoc				

