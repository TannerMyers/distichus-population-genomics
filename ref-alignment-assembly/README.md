
# ***Anolis distichus*** ddRADseq Project

## Assembly of ddRADseq Loci and Calling of Variants -- Mapping Approach

This file outlines workflow of ddRADseq loci assembly via alignment to a reference genome. Individual fastq files are demultiplexed using Stacks' `process_radtags` function before being mapped against the *Anolis carolinensis* reference genome with BWA MEM, then .sam files are converted to .bam and sorted with samtools, assessed for mapping quality with Qualimap, and, finally, SNPs are called using GATK. Final dataset filtering is completed using vcftools following O'Leary et al. (2018).

**Setting up the environment**

I strongly recommend installing packages through a package manager (e.g., Anaconda, Miniconda, or Mamba) 
