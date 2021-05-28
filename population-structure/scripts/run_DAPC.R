####################################################################
# This script performs PCA on genomic SNP data output in STRUCTURE #
# file format output by the populations module of Stacks           #
####################################################################

# Change to your working directory
setwd("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/population-structure/adegenet/")

####################################################################

######################### Load Packages ############################
library(adegenet) # run install.packages() once if not already installed
library(tidyverse)
library(ggplot2)

####################################################################

######################### Load Datasets ############################

# Read structure file in as a dataframe
dat <- read.table("../populations.distichusonly.R0.7.noBahamas.str", 
                  # this file include "Individual" and "Population" labels for first two columns
                  header = FALSE # If header is set to TRUE, include check.names=FALSE
                  )

# Convert structure dataframe to a "genind" object
dat_genind <- df2genind(X = dat,
                        ind.names = dat$V1,
                        pop = dat$V2,
                        sep = '/',
                        NA.char = -9,
                        ploidy = 2,
                        loc.names = dat[1,3:ncol(dat)],
                        strata = NULL,
                        hierarchy = NULL
                        )

# Or, cut out middleman and just directly upload structure file
D <- read.structure("../conStruct/populations.distichusonly.R0.7.noBahamas.str", 
                    # No column names for individual and population columns
                      onerowperind = FALSE, 
                      n.ind=143, 
                      n.loc=1840)
  # Interactive keys: 0, 2, 1, <enter>, 1

####################################################################

################# Evaluate support for values of K #################

# Obtain BIC scores for K values ranging from 1 - 
clusters <- find.clusters(D, 
                     max.n.clust = 10
                    )

# 
table(pop(D), clusters$grp) 
## Note: using locality as population may not be the best way to visualize these
## results as PGDspider encodes each locality in the populations definition file
## as 1, 2, 3, ... for as many sampling localities as you have. If you have a priori
## designations for your populations, you may want to use those to define your populations

# Run with 7 clusters
clusters <- find.clusters(D, clust = NULL, n.clust = 7)

####################################################################

########################### Run DAPC ###############################

results <- dapc(D, pop = clusters$grp, n.pca = 140, scale = FALSE)


