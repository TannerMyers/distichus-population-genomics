setwd("/scratch/tcm0036/distichus-ddRAD/analyses/environmental-association/AnoDist_aln")

rm(list = ls())

# Load R packages
library(tidyverse)
library(vcfR)
library(adegenet)
library(adespatial)
library(sp)
library(rgdal)
library(vegan)

## Load quickMem function from Borcard et al. (2018)
## https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/quickMEM.R
source("~/distichus-ddRAD/scripts/quickMEM.R")


###############################################################
## Anolis distichus subspecies

## Exclude individual specimens that cluster with A. brevirostris
RAD_data <- read_table("~/distichus-ddRAD/info/distichus-subsp-Hispaniola-popmap-cleaned.tsv")

## Load linkage-pruned VCF file with 0% missing data
vcfR <- read.vcfR("/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/eems/AnoDist_aln/distichus/distichus-subsp_variants_no-miss_50KbLD.recode.vcf")
## convert vcfR object to adegenet genind object
dat <- vcfR2genind(vcfR)
	indNames(dat) <- RAD_data$Sample_ID_pop
	pop(dat) <- RAD_data$Locality


## Subset genind object to reflect individuals retained in subsetted information file
    dat <- dat[indNames(dat) %in% RAD_data$Sample_ID_pop]
        print(dat)
        print(dim(dat@tab))

    ## Set population to be sampling locality
    pop(dat) <- RAD_data$Locality

# First, estimate dbMEMs for use as a covariate in RDA
latlong <- SpatialPoints(RAD_data[, c("Longitude", "Latitude")], proj4string = CRS("+proj=longlat"))

distichus_xy <- as.data.frame(spTransform(latlong, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")))

dbmem_1 <- quickMEM(dat@tab, distichus_xy, perm.max = 999)
    save(dbmem_1, file = "rda/dbMEM_distichus-subspecies.RData")

dbmem_2 <- quickMEM(distichus_xy, perm.max = 999)
    save(dbmem_2, file = "rda/dbMEM_2_distichus-subspecies.RData")