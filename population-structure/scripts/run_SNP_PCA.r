####################################################################
# This script performs PCA on genomic SNP data output in STRUCTURE #
# file format output by the populations module of Stacks           #
####################################################################

# Change to your working directory
working_dir <- getwd()
setwd(working_dir)

####################################################################

######################### Load Packages ############################
library(adegenet) # run install.packages() once if not already installed
library(tidyverse)
library(ggplot2)

dat <- read_table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/population-structure/adegenet/populations.structure")

lizards <- df2genind(dat, ploidy = 2)
# Error in df2genind(dat, ploidy = 2) : 
#   Not enough information to convert data: please indicate the separator (sep=...) or the number of characters coding an allele (ncode=...)
