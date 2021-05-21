# ***Anolis distichus*** ddRADseq Project

## Assembly of ddRADseq Loci and Calling of Variants

This file outlines the workflow of demultiplexing and cleaning fastq files, assembling cleaned
and sorted fastq files into loci, and calling variants with the program Stacks. The data are ddRADseq data assembled for
276 individuals of *Anolis distichus* collected across the island of Hispaniola, Hispaniola's satellite
islands Alto Velo and Isla Saona, and the Bahamas and 6 green anoles from the *Anolis carolinensis* species group in Cuba.

Workflow following Rochette & Catchen et al. (2017).

**MODULES TO LOAD:**

- stacks/2.55
- vcftools/v0.1.17
- bwa/0.7.17
- python/3.6.4
- samtools/1.7 
	
	\* samtools and python3 necessary for `stacks-integrate-alignments` 
	
## Demultiplex data with process_radtags:

The first step of stacks is **process_radtags**, which cleans and demultiplexes raw reads kept in fastq files obtained from Novogene.
	
-	ddRADseq libraries were sequenced on 3 Illumina lanes (96 individuals/lane)
- Raw data downloaded from Novogene stored in **Rawdata/** directory
	- **Rawdata/** contains 12 subfolders, each of which contains 6 fastq files
		- P1_1, P1_2, ... P3_4

Run process_radtags by submitting the `demultiplex*_*.sh` scripts corresponding to folders in **Rawdata/**

- Barcode files containing sample IDs are found in **info/** directory
	- Barcodes files for each lane were split into sets of 24, for the number of samples in each subfolder in **Rawdata/**, using the awk one-liner below:
	`awk -v OFS='\t' '{ print $1, $2 }' barcodes_lane3.tsv`
	`> temp && mv temp barcodes_lane3_2.tsv` 
	
- This command within `demultiplex_lane1_1.sh` runs `process_radtags`

		process_radtags 

			-P -1 $filepath1 -2 $filepath2 
			-b $barcodes 
			-o ./ \
			--renz-1 sbfI 
			--renz-2 mspI 
			--inline_null 
			-c -q -r \ 
			&> process_radtags.lane1_1.oe

- Individual sample fastq files are output in the **samples/** directory


Run `processradtags_results.sh` to print standard output and standard error from each process_radtags job, calculate
	summary statistics from reads, and prepare a file for data analysis and visualization in R: `reads_per_sample.tsv`

Then, run R script `proc_radtags_results.R` to plot retained reads per individual


## Choose subset of samples to use for parameter testing:

To identify the combination of parameter values that yields the most loci for population genomic analyses, I used the
wrapper `denovo_map.pl`, on a subset of samples under different values of *M* and *n*.
Following Rochette & Catchen (2017), I created 9 directories in the folder **tests.denovo/** for each value of the 
ustacks parameter, *M*, or the number of mismatches allowed between stacks to merge into a putative locus, named
**stacks.M1**, and so on. I then selected a dozen samples with varying numbers of retained reads following
`process_radtags` and stored them in their own population map file, `popmap.test_samples.tsv`. 

After the stacks pipeline completes, run the module `populations` using the `-R` flag set to 0.8, which requires all loci in the final dataset to occur across 80% of populations. I do this with its own script: `populations.R80.sh` and produce the file `populations.R80.tsv`.

\***********Update after plotting R80 stats \*****************

##Run Stacks for all samples:

Ready to run `denovo_map.pl`, I wrote the shell script, `run_denovomap.sh` that uses a for loop to submit denovo map
for each value of M (1-9) being tested.

	- I used this command to submit the script to the Hopper cluster: 
	`qsub -N denovomaptest -q general -W group_list=jro0014_lab -W x=FLAGS:ADVRES:jro0014_lab -l walltime=48:00:00 run_denovomap.sh` 

- Running `denovo_map.pl` on all samples:
 
	myqsub -N denovomap3all -n -d /scratch/phyletica/distichus/scripts/ --ppn 10 --mem 20gb --time 200:00:00 <<< "./run_denovomap.sh 3"

- Running stacks pipeline manually:

	myqsub -N stacks-pipelineM4n3 -n -d /scratch/phyletica/distichus/scripts/ --ppn 16 --mem 40gb --time 800:00:00 <<< "./run_stacks_pipeline.sh \**insert command line arguments here*\*"
	
	- Memory and time requested should vary. `populations` for the total number of individuals in our dataset needs a lot of memory.
		
	- `Stacks` components included in `run_stacks_pipeline.sh`:

		- `ustacks`
		- `cstacks`
		- `sstacks`
		- `tsv2bam`
		- `gstacks`
		- `bwa mem` to align the consensus catalog loci ('catalog.fa.gz') to the *Anolis carolinensis* reference genome
		- `stacks-integrate-alignments`
		- `populations`
	
	**Making Reference Genome database:**
 
	I downloaded the reference genome for *A. carolinensis* from the Ensembl database with the following line of code:
`rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-101/fasta/anolis_carolinensis/dna/ .`
Check for potential file corruption:

	`sum * > sum.files.txt`

	`diff sum.files.txt CHECKSUMS`

	I used `bwa index` to create a reference genome database:

	Run script `bwa_index.sh` to submit below command as job to Hopper 
	`bwa index -p bwa/anocar $genome_fa &> bwa/bwa_index.oe`
	
	I aligned the consensus sequence of the catalog loci to the *A. carolinensis* reference I had indexed with `bwa mem` in the `run_stacks_pipeline.sh` script and integrated alignment information with the `Stacks` component `stacks-integrate-alignments`.	
	
## Filter individuals:

**Run Stacks on individual populations:**

Following the protocol of Cerca et al. (2021), I ran Stacks on individual groups to generate vcf files that could then be inspected for missing data to identify "bad apples". I generated seven population maps for distichoid anoles, one for each of the following lineages or group of lineages:

* *brevirostris* 
* All Bahamian *subspecies* of *A. distichus*
* *A. distichus* subspecies endemic to the Tiburon penninsula plus *A. d. dominicensis* 3 (Geneva et al. (2015))
* *A. d. dominicensis* 1
* *A. d. dominicensis* 2
* *A. d. dominicensis* 4
* *A. d. ignigularis* 
* *A. d. properus*
* *A. d. sejunctus*
* *A. d. ravitergum* 
* *A. d. favillarum*
* Individuals from localities intermediate between *A. d. dominicensis* 1 and *A. d. dominicensis* 2, *A. d. ravitergum* and *A. d. ignigularis*, and *A. d. ignigularis* and *A. d. properus*. 

Stacks files are output to the **stacks.denovo/population-stacks.denovo/*/** directories	

Run `individual-missing-data-assessment.sh` to loop through directories to run vcftools with the flag `--missing-indv` on the variant call format file "populations.snps.vcf" produced by `populations`, outputting the missing data information to a "bad_apples" file.

I excluded individuals if they exceeded the average percent missing data for their population, unless their percentage of missing data was less than the dataset average of missing data.

**References:**

Cerca, J., M. F. Maurstad, N. C. Rochette, A. G. Rivera‐Colón, N. Rayamajhi, J. M. Catchen, and T. H. Struck. 2021. Removing the bad apples: a simple bioinformatic method to improve loci‐recovery in *de novo* RADseq data for non‐model organisms. Methods Ecol Evol. 

Li, H., and R. Durbin. 2009. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25:1754–1760.

Paris, J. R., J. R. Stevens, and J. M. Catchen. 2017. Lost in parameter space: a road map for stacks. Methods Ecol Evol 8:1360–1373.

Rochette, N. C., A. G. Rivera‐Colón, and J. M. Catchen. 2019. Stacks 2: analytical methods for paired‐end sequencing improve RADseq‐based population genomics. Mol Ecol 28:4737–4754.

Rochette, N. C., and J. M. Catchen. 2017. Deriving genotypes from RAD-seq short-read data using Stacks. Nat Protoc 12:2640–2659.
