#! /bin/bash

# Update once version that works is ID'd
module load stacks 

# Name --popmap variable
popmap=/scratch/phyletica/distichus/info/popmap.test_samples.tsv

# Name --samples variable
samples=/scratch/phyletica/distichus/samples 

# Make an array for the values of M to test
M=( 1 2 3 4 5 6 7 8 9 )

# Name job to be submitted with qsub
jobname=denovomaptest

	# iterate over the array, M
	for M in ${M[@]}; 
			do

		# Name variable for output files 
        out_dir=/scratch/phyletica/distichus/tests.denovo/stacks.M$M 
		# Move into directory for output files
		cd $out_dir

		qsub -N $jobname$M \
                -q general \
                -W group_list=jro0014_lab \
                -W x=FLAGS:ADVRES:jro0014_lab \
                -d $out_dir \
                 -l walltime=48:00:00 << HereDoc

				# run denovo_map for those samples for the value of M being tested
					denovo_map.pl \
							-M $M \
							-n $M \
							-m 3 \
							-o $out_dir \
							--popmap $popmap \
							--samples $samples \
							--paired &> denovomap.oe 	
		
HereDoc			
			done
