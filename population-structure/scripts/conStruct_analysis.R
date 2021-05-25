###############

# Project: Using ddRADseq data to infer patterns of gene flow across the Anolis distichus
# species complex in Hispaniola

# Authors:
# Tanner Myers, Pietro de Mello, Paul Hime, and Rich Glor

# This script estimates population structure while taking spatial  
# information into account to distinguish true sub-structure of   
# populations from isolation-by-distance  

###############

setwd("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/population-structure/conStruct")

# load packages
#library(devtools)
#install_github("gbradburd/conStruct",build_vignettes=TRUE) # Run once
library(conStruct)

myconStructdata <- structure2conStruct(infile = "populations.distichusonly.R0.7.str", 
                                       onerowperind = FALSE, start.loci = 3, 
                                       start.samples = 2, missing.datum = -9, 
                                       outfile = "myconStructdata")

# perform spatial analysis with conStruct
spRun <- conStruct(spatial = TRUE,
                  K =, # update with value of K
                  freqs = myconStructdata,
                  geoDist = ,
                  coords = coords,
                  prefix = "spK") # update with value of K

# perform non-spatial analysis
nspRun <- conStruct(spatial = FALSE,
                    freqs = myconStructdata,
                    geoDist = NULL,
                    coords = coords,
                    prefix = "nspK") # update with value of K
