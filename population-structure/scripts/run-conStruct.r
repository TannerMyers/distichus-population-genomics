setwd("/mmfs1/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/conStruct")

library(tidyverse)
library(fields)
library(conStruct)
library(doParallel)
library(foreach)

#########
# Generate input files for conStruct 
#########

# (1) Load allele frequencies file
    structure2conStruct(infile = "distichus-spgrp-North-Paleo-island-conStruct.str",
        onerowperind = TRUE,
        start.loci = 3,
        missing.datum = 0, # This data matrix lacks any missing data so I just chose 0
        outfile = "distichus-spgrp-North-Paleo-island-conStuct-input")
load("distichus-spgrp-North-Paleo-island-conStuct-input.RData")

# (2) Load coordinates for specimen sampling localities
RAD_data <- read_table("/home/tcm0036/distichus-ddRAD/info/distichus-popmap-cleaned-master.tsv", col_names = TRUE)
RAD_data2 <- RAD_data %>% filter(Sample != "4481" & Sample != "10085" & Sample != "10086") %>%
    filter(!(Taxon %in% c("websteri", "marron", "caudalis", "altavelensis"))) %>%
    filter(Island %in% "Hispaniola")  %>%
    filter(!grepl("_rep", Sample))

snmfK9 <- read_table("../lea/sNMF_K9_ancestry_coefficients_cluster.tsv")
 
RAD_data_north <- RAD_data2[RAD_data2$Sample_ID_pop %in% snmfK9$V1[snmfK9$cluster != 4 & snmfK9$cluster != 6 & snmfK9$cluster != 8],]
latlong <- as.matrix(RAD_data_north[,c('Longitude', 'Latitude')])

# (3) Generate geographic distance matrix from specimen sampling localities

## As per https://cran.r-project.org/web/packages/conStruct/vignettes/format-data.html#geographic-distance-matrix
## I will be using the `fields::rdist.earth` function to estimate pairwise great-circle distances
## between sampling coordinates
geoDist <- fields::rdist.earth(x1 = latlong, x2 =latlong)

#########
# run conStuct analyses for Ks 1-8
#########

dir.create("/mmfs1/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/conStruct/north-island-outputs")

#cl <- makeCluster(10, type = "FORK")
#registerDoParallel(cl)

myRun <- x.validation(train.prop = 0.9, 
    n.reps = 10,
    K = 1:8,
    freqs = freqs,
    coords = latlong,
    geoDist = geoDist,
    n.iter = 1e4,
    make.figs = TRUE,
    save.files = TRUE,
    parallel = TRUE,
    n.nodes = 10, 
    prefix = "north-island-outputs/distichus-North-paleo-island-subsp_CV")

#stopCluster(cl)
