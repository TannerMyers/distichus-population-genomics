setwd("/scratch/tcm0036/distichus-ddRAD/analyses/environmental-association/mmrr")

rm(list = ls())

# Load R packages
library(sp)
library(raster)
library(rgeos)
library(rgdal)
library(RStoolbox)
library(adegenet)
library(vcfR)
library(hierfstat)
library(StAMPP)
library(BEDASSLE)
library(ecodist) # MRM function may be able to replace MMRR
library(gdm)
library(vegan)
library(tidyverse)

# Load MMRR function available from https://datadryad.org/stash/dataset/doi:10.5061%2Fdryad.kt71r
source("~/distichus-ddRAD/scripts/MMRR.R")

# ddRADseq individual information
RAD_data <- read_table("/home/tcm0036/distichus-ddRAD/info/distichus-popmap-cleaned-master.tsv", col_names = TRUE)
RAD_data2 <- RAD_data %>% filter(Sample != "4481" & Sample != "10085" & Sample != "10086") %>%
    filter(!(Taxon %in% c("websteri", "marron", "caudalis", "altavelensis"))) %>%
    filter(Island %in% "Hispaniola")  %>%
    filter(!grepl("_rep", Sample)) # Temporary, eventually I need to deal with the technical replicates in the dataset

# ddRADseq sequence data
vcfR <- read.vcfR("../../population-structure/distichus-spgrp-no-missing-LD-pruned-informative-names.vcf")

# convert vcfR object to genind object
    dat <- vcfR2genind(vcfR)
    dat2 <- dat[indNames(dat) %in% RAD_data2$Sample_ID_pop]
    pop(dat2) <- RAD_data2$Locality

# Estimate pairwise Fst for genetic differentiation matrix
data <- genind2hierfstat(dat = dat2, pop = dat2@pop)

gen_diff <- pairwise.WCfst(dat = data, diploid = TRUE)
	save(gen_diff, file = "Fst-matrix.RData")
