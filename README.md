# ***Anolis distichus*** ddRADseq Project

This file outlines the workflow of demultiplexing and cleaning fastq files, assembling cleaned
and sorted fastq files into loci, and calling variants with the program Stacks. The data are ddRADseq data assembled for
276 individuals of *Anolis distichus* collected across the island of Hispaniola, Hispaniola's satellite
islands Alto Velo and Isla Saona, and the Bahamas and 6 green anoles from the *Anolis carolinensis* species group in Cuba.

Workflow following Rochette & Catchen et al. (2017).

**MODULES TO LOAD:**

- stacks/2.5
- bwa/0.7.17*
 	* bwa necessary only for assembling reference genome database

**process_radtags:**

The first step of stacks is **process_radtags**, which cleans and demultiplexes raw reads kept in fastq files obtained from Novogene.
	
-	ddRADseq libraries were sequenced on 3 Illumina lanes (96 individuals/lane)
- Raw data downloaded from Novogene stored in **Rawdata/** directory
	- **Rawdata/** contains 12 subfolders, each of which contains 6 fastq files
		- P1_1, P1_2, ... P3_4

Run process_radtags by submitting the `demultiplex*_*.sh` in the **cleaned** directories corresponding to folders in **Rawdata/**

- Barcode files containing sample IDs are found in **info/** directory
	- Barcodes files for each lane were split into sets of 24, for the number of samples in each subfolder in **Rawdata/**, using the awk one-liner below:
	`awk -v OFS='\t' '{ print $1, $2 }' barcodes_lane3.tsv`
	`> temp && mv temp barcodes_lane3_2.tsv` 
	
- This command within `demultiplex_lane1_1.sh` runs **process_radtags**

`process_radtags -P -1 $filepath1 -2 $filepath2 -b $barcodes -o ./ \`

`--renz-1 sbfI --renz-2 mspI --inline_null -c -q -r \` 
`&> process_radtags.lane1_1.oe`

- Individual sample fastq files are output in the appropriate cleaned directory


Run `processradtags_results.sh` to print standard output and standard error from each process_radtags job, calculate
	summary statistics from reads, and prepare a file for data analysis and visualization in R: `reads_per_sample.tsv`

Then, run R script `proc_radtags_results.R` to plot retained reads per individual


**Choose subset of samples to use for parameter testing:**

To identify the combination of parameter values that yields the most loci for population genomic analyses, I used the
wrapper Perl script, `denovo_map.pl`, on a subset of samples under different values of *M*, *m*, and *n*.
Following Rochette & Catchen (2017), I created 9 directories in the folder **tests.denovo/** for each value of the 
ustacks parameter, *M*, or the number of mismatches allowed between stacks to merge into a putative locus, named
**stacks.M1**, **stacks.M2**, and so on. I then selected six samples with varying numbers of retained reads following
**process_radtags** and stored them in their own population map file, `popmap.test_samples.tsv`. 

Ready to run `denovo_map.pl`, I wrote the shell script, `run_denovomap.sh` that uses a for loop to submit denovo map
for each value of M (1-9) being tested.

	- I used this command to submit the script to the Hopper cluster: 
	`qsub -N denovomaptest -q general -W group_list=jro0014_lab -W x=FLAGS:ADVRES:jro0014_lab -l walltime=48:00:00 run_denovomap.sh` 

- Running `denovo_map.pl` on all samples:
 
	myqsub -N denovomap3all -n -d /scratch/phyletica/distichus/scripts/ --ppn 10 --mem 20gb --time 200:00:00 <<< "./run_denovomap.sh 3"



**Making Reference Genome database:**
 
I downloaded the reference genome for *A. carolinensis* from the Ensembl database with the following line of code:
`rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-101/fasta/anolis_carolinensis/dna/ .`
Check for potential file corruption:

`sum * > sum.files.txt`

`diff sum.files.txt CHECKSUMS`

Following Rochette & Catchen (2017), I used `bwa index` to create a reference genome database:

Run script `bwa_index.sh` to submit below command as job to Hopper 
	`bwa index -p bwa/anocar $genome_fa &> bwa/bwa_index.oe`
	
Once 	
**References:**

Li, H., and R. Durbin. 2009. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25:1754–1760.

Rochette, N. C., and J. M. Catchen. 2017. Deriving genotypes from RAD-seq short-read data using Stacks. Nat Protoc 12:2640–2659..
