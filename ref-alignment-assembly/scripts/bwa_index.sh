#!/usr/bin/env bash
#SBATCH --job-name=bwa-index
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tcm0036@auburn.edu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task 8
#SBATCH --mem 32G
#SBATCH --partition jro0014_amd

# Assign variables
## Download genome from Ensembl FTP site beforehand
# genome_fa=/home/tcm0036/distichus-ddRAD/genome/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa
# bwa_db=/home/tcm0036/distichus-ddRAD/genome/bwa/anocar
genome_fa=/home/tcm0036/distichus-ddRAD/genome/Anolis_distichus/AnoDis1.0.fasta.gz
bwa_db=/home/tcm0036/distichus-ddRAD/genome/Anolis_distichus/bwa/anodis

# Load bwa module
module load bwa

# Create index for reference genome
# bwa index -p $bwa_db $genome_fa ## Not sure why this doesn't work
bwa index $genome_fa