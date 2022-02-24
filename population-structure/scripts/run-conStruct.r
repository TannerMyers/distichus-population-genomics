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
    #structure2conStruct(infile = "distichus-spgroup-noheader-conStruct.str",
    #    onerowperind = TRUE,
    #    start.loci = 3,
    #    missing.datum = 0, # This data matrix lacks any missing data so I just chose 0
    #    outfile = "distichus-spgroup-conStuct-input")
load("distichus-spgroup-conStuct-input.RData")

# (2) Load coordinates for specimen sampling localities
ind_data <- read_table("/home/tcm0036/distichus-ddRAD/info/all-distichoids-Hispaniola-only-cleaned-popmap.tsv", col_names = TRUE)
latlong <- as.matrix(ind_data[,c('Longitude', 'Latitude')])

# (3) Generate geographic distance matrix from specimen sampling localities

## As per https://cran.r-project.org/web/packages/conStruct/vignettes/format-data.html#geographic-distance-matrix
## I will be using the `fields::rdist.earth` function to estimate pairwise great-circle distances
## between sampling coordinates
geoDist <- fields::rdist.earth(x1 = latlong, x2 =latlong)

#########
# run conStuct analyses for Ks 1-10
#########

dir.create("/mmfs1/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/conStruct/outputs")

#cl <- makeCluster(10, type = "FORK")
#registerDoParallel(cl)

myRun <- x.validation(train.prop = 0.9, 
    n.reps = 10,
    K = 1:10,
    freqs = freqs,
    coords = latlong,
    geoDist = geoDist,
    n.iter = 1e4,
    make.figs = TRUE,
    save.files = TRUE,
    parallel = TRUE,
    n.nodes = 10, 
    prefix = "outputs/distichus-spgroup-CV")

#stopCluster(cl)
