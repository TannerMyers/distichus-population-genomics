#! /usr/bin/env bash

################################################################

# This script uses GATK to call variants and generate a multi-sample 
# VCF file.

#Note: Prior to running this script, users should create a 
# conda environment with the packages involved here if they 
# are not using their HPC's modules

################################################################

# Name variables 
BAM_DIR=/scratch/phyletica/distichus/alignments/bam-files/

VCF_DIR=/scratch/phyletica/distichus/alignments/GATK/vcf/


# Assign variable for reference genome fasta
GENOME=/scratch/phyletica/distichus/genome/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa

# Designate sample ID variables

    # Population map file containing individual IDs that for loop will iterate over.
    POP_MAP=/scratch/phyletica/distichus/info/popmap.tsv
    # Isolate individual IDs in population map to use for looping variable
    SAMPLES=`awk 'BEGIN {OFS = FS} {print $1}' $POP_MAP`

################################################################

# Run HaplotypeCaller
for sample in $SAMPLES; 
    do
        # For each sample, submit a job to align the reads to the indexed genome with BWA MEM, convert to .bam and sort with samtools,
        # and perform mapping QC with Qualimap.
        qsub -N $sample.GATK -d /scratch/phyletica/distichus/scripts -q gen28 -W group_list=jro0014_lab -W x=FLAGS:ADVRES:jro0014_s28 -l nodes=1:ppn=8,mem=80gb,walltime=12:00:00 <<<" 
        # Load conda environment with bwa, samtools, and qualimap
        source ~/mambaforge/etc/profile.d/conda.sh # Or path to where your conda is
        conda activate genomics_env

        gatk --java-options "-Xmx4g" HaplotypeCaller \
                --reference $GENOME \
                --input $BAM_DIR$sample.sorted.bam \
                --output $VCF_DIR$sample.g.vcf.gz \
                --emitRefConfidence GVCF
        "
        done 
