#! /usr/bin/env bash

#source ~/mambaforge/etc/profile.d/conda.sh
#conda activate eems_env


for i in 1 2 3;
	do

		sbatch --job-name=EEMS-chain$i --partition=jro0014_amd --cpus-per-task=1 --mem=60G --time=30-24:00:00 --wrap="
			source ~/mambaforge/etc/profile.d/conda.sh
			conda activate eems_env
									
				runeems_snps --params /home/tcm0036/distichus-ddRAD/scripts/params-chain${i}.ini
			"

	done
