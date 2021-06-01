###############

# Project: Using ddRADseq data to infer patterns of gene flow across the Anolis distichus
# species complex in Hispaniola

# Authors:
# Tanner Myers, Pietro de Mello, Paul Hime, and Rich Glor

# This script estimates population structure while taking spatial  
# information into account to distinguish true sub-structure of   
# populations from isolation-by-distance  

###############

# Change to your working directory
setwd("/gpfs01/home/tcm0036/distichus/population-structure/conStruct")

# Load packages 


## Include the "lib.loc" argument is necessary to tell R on the cluster where to look for the installed package
library(optparse)
library(gtools, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(doParallel, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(Rcpp, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(RcppParallel, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(RcppEigen, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(rstan, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library") 
# the lib.loc argument is necessary to tell R on the cluster where to look for the installed package
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
#library(devtools)
#install_github("gbradburd/conStruct",build_vignettes=TRUE) # Run once
library(conStruct, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")


# Load input files 

# Load allele frequency data
#myconStructdata <- structure2conStruct(infile = "../populations.distichusonly.R0.7.noBahamas.str", 
#                                       onerowperind = FALSE, start.loci = 3, 
#                                       start.samples = 2, missing.datum = -9, 
#                                       outfile = "myconStructdata")
## Instead of using this function from conStruct, process the STRUCTURE file you obtain
## from Stacks using the "structure-to-conStruct.R" script employing Emily Puckett's code


# Use processed structure file obtained by running "structure-to-conStruct.R"
freq <- read.table("populations.distichusonly.R0.7.noBahamas.processed.str")
freq <- as.matrix(freq)

# Load geographic distance matrix generated with topoDistance package
load(file = "ddRAD_topoDistance_geog_mat.rda")

# Load lat. and long. for sampling localities
latlong <- read.csv("ddRAD_cleaned_Hispaniolan_distichus_localities.csv", # Update path here
                   header=TRUE)
latlong <- as.matrix(latlong[, c('Longitude', 'Latitude')])
class(latlong)

# Loop over values of K, performing both spatial and non-spatial runs
          

#for(k in 1:10){

# perform spatial analysis with conStruct
spRun.K <- conStruct(spatial = TRUE,
                K = as.numeric(k), # update with value of K
                freqs = freq,
                geoDist = tDist.mat,
                coords = latlong,
		        n.chains = 8,
		        n.iter = 10000,
	            save.files = TRUE,
		        make.figs = FALSE,
                prefix = paste0("spK",k)) # update with value of K

# perform non-spatial analysis
nspRun.K <- conStruct(spatial = FALSE,
                K = as.numeric(k),
		        freqs = freq,
                geoDist = NULL,
                coords = latlong,
		        n.chains = 8,
		        n.iter = 10000,
		        save.files = TRUE,
		        make.figs = FALSE,
                prefix = paste0("nspK",k)) # update with value of K

}

