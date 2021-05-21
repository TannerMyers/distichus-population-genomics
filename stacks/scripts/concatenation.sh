#! /bin/bash

# "Cleaned" is the directory where the subfolders containing all fastq files
# output by `process_radtags` are located
file_dir=/scratch/phyletica/distichus/cleaned
	cd $file_dir

# Make array based on how lanes were divided by sequencer e.g., P1_1
dirs=(1_1 1_2 1_3 1_4 2_1 2_2 2_3 2_4 3_1 3_2 3_3 3_4)

# Designate directory to contain concatenated sample fastqs
out_dir=/scratch/phyletica/distichus/cleaned/samples

for dir in ${dirs[@]}; 
	do 
	# Move into the relevant folder with processed fastq files	
		cd cleaned$dir
	
	# Specify the barcodes file for that subfolder
		barcodes=/scratch/phyletica/distichus/info/barcodes_lane$dir.tsv
	
	# Make an array out of the samples included in the barcodes file 
		samples=($(awk 'FS="\t" {print $2}' $barcodes))
		
			# Iterate through the just made array to concatenate all of the
			# fastq files for a sample into a single fastq file
			for sample in "${samples[@]}"; do
				cp $sample.1.fq.gz $out_dir/$sample.fq.gz
                       		cat $sample.2.fq.gz >> $out_dir/$sample.fq.gz
                        	cat $sample.rem.1.fq.gz >> $out_dir/$sample.fq.gz
                       		cat $sample.rem.2.fq.gz >> $out_dir/$sample.fq.gz
			done

	# Move back to the main "cleaned" directory
		cd $file_dir
	done
