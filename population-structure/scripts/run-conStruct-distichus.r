setwd("/mmfs1/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/conStruct/AnoDist_aln/distichus-subsp")

library(tidyverse)
#library(fields)
library(sp)
library(raster)
library(conStruct)
library(doParallel)
library(foreach)

print(Sys.time())

#########
# Generate input files for conStruct
#########

# (1) Load allele frequencies file
    # structure2conStruct(infile = "distichus-subsp.recode.strct_in",
    #     onerowperind = TRUE,
    #     start.loci = 3,
    #     start.samples = 3,
    #     missing.datum = 0, # This data matrix lacks any missing data so I just chose 0
    #     outfile = "distichus-conStruct-input")
load("distichus-conStruct-input.RData")

# (2) Load coordinates for specimen sampling localities
RAD_data <- read_table("~/distichus-ddRAD/info/distichus-subsp-Hispaniola-popmap-cleaned.tsv", col_names = TRUE)
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

#myRun <- x.validation(train.prop = 0.8,
#    n.reps = 10,
#    K = 1:10,
#    freqs = freqs,
#    coords = latlong,
#    geoDist = geog_dist,
#    n.iter = 1e4,
#    make.figs = TRUE,
#    save.files = TRUE,
#    parallel = TRUE,
#    n.nodes = 10,
#    prefix = "outputs/distichus-subsp-")
#warnings()

## Unfortunately, as of 2/27/2023, this job ended before completing
## so I will attempt to load the data.partitions object created in the
## first attempt so I don't have to start from scratch
# conStruct:::x.validation.rep
# function (rep.no, K, data.partition, geoDist, coords, prefix, n.iter, make.figs = FALSE, save.files = FALSE, ...) 

load("outputs/distichus-subsp-.xval.data.partitions.Robj")
conStruct:::x.validation.rep(data.partition = data.partitions[[1]],
        prefix = "outputs/distichus-subsp-",
        save.files = FALSE,
        make.figs = FALSE,
        n.iter = 1e4,
        K = 1:10,
        geoDist = geog_dist,
        coords = latlong)

## Following Gideon's response to a GitHub issue
## see here: https://github.com/gbradburd/conStruct/issues/40
## I will run the following code to standardize cross-validation
## results and write out the output.
#names(x.val) <- paste0("rep_", 1:n.reps)
#x.val <- lapply(x.val, conStruct:::standardize.xvals)
#save(x.val, file = paste0("outputs/distichus-subsp", ".xval.results.Robj")
#conStruct:::write.xvals(x.val, "outputs/distichus-subsp")

#stopCluster(cl)
