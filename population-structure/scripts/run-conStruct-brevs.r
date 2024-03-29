setwd("/mmfs1/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/conStruct/AnoDist_aln/brevirostris")

library(tidyverse)
#library(fields)
library(sp)
library(raster)
library(conStruct)
library(doParallel)
library(foreach)

print("The time is")
Sys.time()

#########
# Generate input files for conStruct
#########

# (1) Load allele frequencies file
# structure2conStruct(infile = "brevirostris.recode.strct_in",
#     onerowperind = TRUE,
#     start.loci = 3,
#     start.samples = 3,
#     missing.datum = 0, # This data matrix lacks any missing data so I just chose 0
#     outfile = "brevirostris-conStruct-input")
load("brevirostris-conStruct-input.RData")

# (2) Load coordinates for specimen sampling localities
RAD_data <- read_table("~/distichus-ddRAD/info/brevirostris-popmap-cleaned.tsv", col_names = TRUE)
latlong <- as.matrix(RAD_data[,c('Longitude', 'Latitude')])

# (3) Generate geographic distance matrix from specimen sampling localities

## As per https://cran.r-project.org/web/packages/conStruct/vignettes/format-data.html#geographic-distance-matrix
## I will be using the `fields::rdist.earth` function to estimate pairwise great-circle distances
## between sampling coordinates
# geoDist <- fields::rdist.earth(x1 = latlong, x2 =latlong) # Not symmetric for some reason
geog_dist <- sp::spDists(x = latlong, y = latlong, longlat = TRUE)

#########
# run conStuct analyses for Ks 1-8
#########

#cl <- makeCluster(10, type = "FORK")
#registerDoParallel(cl)

myRun <- x.validation(train.prop = 0.8,
    n.reps = 10,
    K = 1:4,
    freqs = freqs,
    coords = latlong,
    geoDist = geog_dist,
    n.iter = 5e3,
    make.figs = TRUE,
    save.files = TRUE,
    prefix = "outputs/second-attempt/brevirostris-"
    # parallel = TRUE,
    # n.nodes = 10,
)

#stopCluster(cl)