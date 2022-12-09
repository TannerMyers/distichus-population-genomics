#!/usr/bin/env bash
#SBATCH --job-name variant-calling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tcm0036@auburn.edu
#SBATCH --time=3-24:00:00
#SBATCH --cpus-per-task 8 
#SBATCH --mem 80G
#SBATCH --partition jro0014_amd

# Must be bgzipped if compressed
genome=/home/tcm0036/distichus-ddRAD/genome/Anolis_distichus/AnoDis1.0.fasta
popmap=/home/tcm0036/distichus-ddRAD/info/popmap-bcftools.tsv
outputfile=/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/bcf/variants.bcf
tmp=mapped-files.tmp.list

ls /scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/bam/*sorted.bam > $tmp 

module load bcftools

bcftools mpileup -Ou -a AD --fasta-ref $genome --bam-list $tmp | \
bcftools call --variants-only --threads $SLURM_CPUS_PER_TASK --multiallelic-caller --group-samples $popmap -Ob -o $outputfile

rm $tmp