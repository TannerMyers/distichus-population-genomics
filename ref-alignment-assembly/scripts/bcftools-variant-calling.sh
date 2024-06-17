#!/usr/bin/env bash
#SBATCH --job-name variant-calling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tcm0036@auburn.edu
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task 16
#SBATCH --mem 30G
#SBATCH --partition bac0071_amd

# Load conda environment containing R
source ~/mambaforge/etc/profile.d/conda.sh # Or path to where your conda is
conda activate genomics_env

# Must be bgzipped if compressed
genome=/home/tcm0036/distichus-ddRAD/genome/Anolis_distichus/AnoDis1.0.fasta
popmap=/home/tcm0036/distichus-ddRAD/info/distichus-popmap-bcftools.tsv
outputfile=/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/bcf/distichus_variants.bcf 
tmp=mapped-files.list.tmp

awk '{print $1}' $popmap > $tmp

bcftools mpileup -Ou --annotate AD,DP --fasta-ref $genome --bam-list $tmp | \
bcftools call --variants-only --threads $SLURM_CPUS_PER_TASK --multiallelic-caller --group-samples $popmap -Ob -o $outputfile

#rm $tmp