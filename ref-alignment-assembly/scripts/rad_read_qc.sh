#!/usr/bin/env bash
#SBATCH --job-name=rad_qc
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=tcm0036@auburn.edu
#SBATCH --time=24:00:00
#SBATCH --mem 30G
#SBATCH --partition bac0071_amd

# Name variables


popmap=/home/tcm0036/distichus-ddRAD/info/distichus-popmap-master.tsv

samples=`awk 'BEGIN {OFS = FS} {print $1}' $popmap | tail -n283`

in_dir=/home/tcm0036/distichus-ddRAD/samples
out_dir=/home/tcm0036/distichus-ddRAD/samples/trimmed

[[ -d $out_dir ]] || mkdir -p $out_dir 

# Load conda environment with fastp
    source ~/mambaforge/etc/profile.d/conda.sh # Or path to where your conda is
    conda activate genomics_env

for ID in $samples;
	do
		
		infq1=$in_dir/$ID.1.fq.gz
		infq2=$in_dir/$ID.2.fq.gz
		outfq1=$out_dir/$ID.1.fq.gz
		outfq2=$out_dir/$ID.2.fq.gz

		fastp --in1 $infq1 --in2 $infq2 --out1 $outfq1 --out2 $outfq2 -l 100 -h distichus-ddRADseq_qc.html &> radqc.log
	
	done
