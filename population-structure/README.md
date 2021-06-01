# ***Anolis distichus*** ddRADseq Project

## Population Structure Analyses

### Dealing with *Stacks* outputs

For population structure analyses, you should really only need the following output files from Stacks (and can get by with as few as one): a .vcf (--vcf argument of populations), the plink format .ped and .map files (--plink argument of populations), or the STRUCTURE output file (--structure argument of populations). However, these files aren't ready to be plugged into an analysis and you'll need to do some editing/converting to get them ready.

- First, some programs may not tolerate the first line of Stacks output files that contains a timestamp of when they were generated. I have copied the original versions of these files and removed this line for input in programs like plink, Admixture, etc.

- To prepare the `plink` output of `populations` for Admixture (same can be done with the `vcf` output by providing plink with the `--vcf` argument.

	- First convert the .ped & .map files to .bed file format (binary plink). The `--allow-extra-chr` argument is necessary when dealing with non-chromosome level data (i.e., data with locations that are on scaffolds)

			plink --file populations.plink --allow-extra-chr
		 	--no-sex --make-bed 
		 	--out 
		 	plink/distichus-only/
		 	populations.cleaned.R0.7.nobrevirostiris
		
	- Then, convert back to .ped

			plink --bfile 
			populations.cleaned.R0.7.nobrevirostiris
			--allow-extra-chr --recode12 
			--out plink/distichus-only/
			populations.cleaned.R0.7.nobrevirostiris
			
			
	
### PCA

To get PCA results, all you need to do is add the `--pca` flag to a plink command, which will result in two output files with the extensions .eigenvals and .eigenvec containing eigenvalues and eigenvectors respectively. These can then be loaded into R and plotted. I did so with the script: `plink-PCA-results.R`.

### Admixture

Admixture v.1.3.0 is parametric method for inferring population structure using a fast maximum-likelihood algorithm. I used Admixture to assign individuals in my dataset to population clusters under a range of K values and then estimated support using its cross validation procedure.

- First, if you don't have it already available, install Admixture using conda:

		conda install -c bioconda admixture

- Then, run Admixture by supplying `admixture` with your .ped and .map files or .bed file and desired value of K. I used the `run_admixture.sh` script to perform many Admixture runs by looping over different K values.

- The optimal value of K can be found in the log*.out files

- Admixture results can be plotted several different ways: 
	- `plot_Admixture.r` uses functions from the ggplot2 and tidyverse R packages to make barplots from the .Q files output by Admixture
	- `admixture-proportions-plotting.R` uses the `make.admix.pie.plot` and `make.structure.plot` functions from the conStruct R package to plot Admixture results as piecharts corresponding to sampling localities as well as barplots.


### conStruct