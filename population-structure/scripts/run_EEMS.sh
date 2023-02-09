#! /usr/bin/env bash

for species in distichus brevirostris;
	do
		for i in 200 400;
			do
				for j in 1 2 3;
					do
						sbatch --job-name=${species}-EEMS-nDemes${i}-chain${j} --partition=bac0071_amd --mem=60G --time=7-24:00:00 --wrap="

							source ~/mambaforge/etc/profile.d/conda.sh
							conda activate eems_env

							echo /scratch/tcm0036/distichus-ddRAD/analyses/population-structure/eems/AnoDist_aln/${species}/param-files/params-chain${j}-${i}demes.ini
							runeems_snps --params /scratch/tcm0036/distichus-ddRAD/analyses/population-structure/eems/AnoDist_aln/${species}/param-files/params-chain${j}-${i}demes.ini
							"
					done

			done
	done