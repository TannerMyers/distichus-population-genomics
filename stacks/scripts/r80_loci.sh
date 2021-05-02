#!/usr/bin/env bash 

M=$1

STACKS_DIR=/scratch/phyletica/distichus/tests.denovo/stacks.M$M

OUT_DIR=$STACKS_DIR/populations.r80

LOG=$OUT_DIR/populations.oe

JOB_NAME=populations-r80

	        qsub -N $JOB_NAME$M \
                         -q general \
                         -W group_list=jro0014_lab \
                         -W x=FLAGS:ADVRES:jro0014_lab \
                         -d /scratch/phyletica/distichus/tests.denovo \
                         -l walltime=48:00:00 << HereDoc

			module load stacks
	
			populations \
				-P $STACKS_DIR \
				-O $OUT_DIR \
				-r 0.80 &> $LOG

HereDoc
		
