#! /bin/bash

module load stacks

# tell denovo map the samples you want to test parameters with
popmap=/scratch/phyletica/distichus/info/popmap.test_samples.tsv

# tell denovo map where the fastq files for those samples are
reads_dir=/scratch/phyletica/distichus/cleaned/Test

# make an array for the values of the Stacks parameter M you
# want to test
M=( 1 2 3 4 5 6 7 8 9 )

	# iterate over the array, M
	for M in ${M[@]}; 
			do
		
		# where denovo map results should be deposited  
		out_dir=/scratch/phyletica/distichus/tests.denovo/stacks.M$M
			cd $out_dir
		
		# run denovo_map for those samples for the value of M being tested
		## Code obtained from Rochette & Catchen (2017) Nature Protocols
		denovo_map.pl --samples $reads_dir --popmap $popmap -o $out_dir \
				-M $M -n $M -m 3 --paired &> denovomap.oe 	
		# removed "-b 1" and "-S" flags as it does not appear to be present in Stacks v2.5
		# -b flag supposed to identify the batch ID for the dataset in question
		# -S flag disables recording of SQL data in database
	
			done
