setwd("/scratch/tcm0036/distichus-ddRAD/analyses/environmental-association/")

rm(list = ls())

# Load R packages
library(tidyverse)
library(reshape2)
library(vcfR)
library(adegenet)
library(hierfstat)
library(sp)
library(raster)
library(rgeos)
library(rgdal)
library(adespatial)
library(envirem)
library(RStoolbox)
library(BEDASSLE)
library(ecodist) # MRM function may be able to replace MMRR
library(gdm)
library(vegan)
library(lavaan)
#library(devtools)
#devtools::install_github("cran/simba")
library(simba)


## Load quickMem function from Borcard et al. (2018)
## https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/quickMEM.R
source("~/distichus-ddRAD/scripts/quickMEM.R")

# Load genetic and morphological datasets and sampling info

# Load ddRADseq individual information
RAD_data <- read_table("/home/tcm0036/distichus-ddRAD/info/distichus-popmap-cleaned-master.tsv", col_names = TRUE)
## Subset dataset to include only Hispaniolan members of the distichus species group
RAD_data2 <- RAD_data %>% filter(Sample != "4481" & Sample != "10085" & Sample != "10086") %>%
    filter(!(Taxon %in% c("websteri", "marron", "caudalis", "altavelensis"))) %>%
    filter(!grepl("_rep", Sample)) %>% # Temporary, eventually I need to deal with the technical replicates in the dataset
    filter(Island %in% "Hispaniola")


snmfK9 <- read_table("/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/lea/sNMF_K9_ancestry_coefficients_cluster-specimen-cleaned.tsv")

# Load linkage-pruned VCF file with 0% missing data
vcfR <- read.vcfR("../population-structure/distichus-spgrp-no-missing-LD-pruned-informative-names.vcf")

#############################################################
# Anolis distichus + A. brevirostris 

# convert vcfR object to adegenet genind object
    dat <- vcfR2genind(vcfR)
    ## Subset genind object to reflect individuals retained in subsetted information file
#    dat2 <- dat[indNames(dat) %in% RAD_data2$Sample_ID_pop]
    ## Set population to be sampling locality
#    pop(dat2) <- RAD_data2$Locality

##############################################################
# First, estimate dbMEMs for use as a covariate in RDA
#latlong <- SpatialPoints(RAD_data2[, c("Longitude", "Latitude")], proj4string = CRS("+proj=longlat"))

#distichoid_xy <- as.data.frame(spTransform(latlong, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")))

#dbmem.1 <- quickMEM(dat2@tab, distichoid_xy, perm.max = 999)
#    save(dbmem.1, file = "rda/dbMEM_1.RData")

#dbmem.2 <- quickMEM(distichoid_xy, perm.max = 999)
#    save(dbmem.2, file = "rda/dbMEM_2.RData")


###############################################################
# Anolis distichus subspecies

## Exclude individual specimens that cluster with A. brevirostris
RAD_data3 <- RAD_data2[RAD_data2$Sample_ID_pop %in% snmfK9$V1[snmfK9$cluster != 6 & snmfK9$cluster != 8], ]

## Subset genind object to reflect individuals retained in subsetted information file
    dat3 <- dat[indNames(dat) %in% RAD_data3$Sample_ID_pop]
        print(dat3)
        print(dim(dat3@tab))

    ## Set population to be sampling locality
    pop(dat3) <- RAD_data3$Locality

# First, estimate dbMEMs for use as a covariate in RDA
latlong <- SpatialPoints(RAD_data3[, c("Longitude", "Latitude")], proj4string = CRS("+proj=longlat"))

distichus_xy <- as.data.frame(spTransform(latlong, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")))

dbmem.1 <- quickMEM(dat3@tab, distichus_xy, perm.max = 999)
    save(dbmem.1, file = "rda/dbMEM_distichus-subspecies.RData")

dbmem.2 <- quickMEM(distichus_xy, perm.max = 999)
    save(dbmem.2, file = "rda/dbMEM_2_distichus-subspecies.RData")