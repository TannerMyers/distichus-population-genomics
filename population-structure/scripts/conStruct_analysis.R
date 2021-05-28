###############

# Project: Using ddRADseq data to infer patterns of gene flow across the Anolis distichus
# species complex in Hispaniola

# Authors:
# Tanner Myers, Pietro de Mello, Paul Hime, and Rich Glor

# This script estimates population structure while taking spatial  
# information into account to distinguish true sub-structure of   
# populations from isolation-by-distance  

###############

#setwd("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/population-structure/conStruct")

# Load packages
# Include the "lib.loc" argument is necessary to tell R on the cluster where to look for the installed package
library(optparse, lib.loc="miniconda3/envs/R_env/lib/R/library/")
library(doParallel, lib.loc="miniconda3/envs/R_env/lib/R/library/")
library(Rcpp, lib.loc="miniconda3/envs/R_env/lib/R/library/")
library(RcppParallel, lib.loc="miniconda3/envs/R_env/lib/R/library/")
library(RcppEigen, lib.loc="miniconda3/envs/R_env/lib/R/library/")
library(rstan, lib.loc="miniconda3/envs/R_env/lib/R/library/") 
# the lib.loc argument is necessary to tell R on the cluster where to look for the installed package
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
#library(devtools)
#install_github("gbradburd/conStruct",build_vignettes=TRUE) # Run once
library(conStruct, lib.loc="miniconda3/envs/R_env/lib/R/library/")

# Load allele frequency data
#myconStructdata <- structure2conStruct(infile = "../populations.distichusonly.R0.7.noBahamas.str", 
#                                       onerowperind = FALSE, start.loci = 3, 
#                                       start.samples = 2, missing.datum = -9, 
#                                       outfile = "myconStructdata")

# Use processed structure file obtained by running ".R"
freq <- read.table("populations.distichusonly.R0.7.noBahamas.processed.str")
freq <- as.matrix(freq)

# Load geographic distance matrix generated with topoDistance package
load(file = "distance-matrix/ddRAD_topoDistance_geog_mat.rda") # Update path here

# Load lat. and long. for sampling localities
latlong <- read.csv("ddRAD_cleaned_Hispaniolan_distichus_localities.csv", # Update path here
                   header=TRUE)
latlong <- as.matrix(latlong[, c('Longitude', 'Latitude')])
class(latlong)

# perform spatial analysis with conStruct
spRun <- conStruct(spatial = TRUE,
                  K =7, # update with value of K
                  freqs = freq,
                  geoDist = tDist.mat,
                  coords = latlong,
                  prefix = "spK7") # update with value of K

# perform non-spatial analysis
nspRun <- conStruct(spatial = FALSE,
                    freqs = freq,
                    geoDist = NULL,
                    coords = coords,
                    prefix = "nspK") # update with value of K
