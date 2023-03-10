setwd("/mmfs1/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/conStruct/AnoDist_aln/distichus-subsp")

library(tidyverse)
library(sp)
library(raster)
library(conStruct)
library(doParallel)
library(foreach)

print(Sys.time())

#########
## Generate input files for conStruct
#########

## (1) Load allele frequencies file
    # structure2conStruct(infile = "distichus-subsp.recode.strct_in",
    #     onerowperind = TRUE,
    #     start.loci = 3,
    #     start.samples = 3,
    #     missing.datum = 0, # This data matrix lacks any missing data so I just chose 0
    #     outfile = "distichus-conStruct-input")
# load("distichus-conStruct-input.RData")

structure2conStruct(infile = "north-island/north-island-subsp.recode.strct_in", 
	 onerowperind = TRUE,  
	 start.loci = 3, 
	 start.samples = 3, 
 	 missing.datum = 0,  
	 outfile = "north-island/north-subsp-input")
load("north-island/north-subsp-input.RData")

structure2conStruct(infile = "south-island/south-island-subsp.recode.strct_in",
	 onerowperind = TRUE, 
	 start.loci = 3, 
	 start.samples = 3, 
	 missing.datum = 0, 
	 outfile = "south-island/south-subsp-input")
load("south-island/south-subsp-input.RData")

# (2) Load coordinates for specimen sampling localities

# RAD_data <- read_table("~/distichus-ddRAD/info/distichus-subsp-Hispaniola-popmap-cleaned.tsv", col_names = TRUE)
# latlong <- as.matrix(RAD_data[,c('Longitude', 'Latitude')])

## North paleo-island pop map
RAD_data <- read_table("~/distichus-ddRAD/info/distichus-north-island-popmap-master.tsv", col_names = TRUE)
latlong <- as.matrix(RAD_data[,c('Longitude', 'Latitude')])

## South paleo-island pop map
RAD_data <- read_table("~/distichus-ddRAD/info/distichus-south-island-popmap-master.tsv", col_names = TRUE)
latlong <- as.matrix(RAD_data[,c('Longitude', 'Latitude')])

## (3) Generate geographic distance matrix from specimen sampling localities

#geoDist <- fields::rdist.earth(x1 = latlong, x2 =latlong) # What Gideon Bradburd uses in vignette, but not symmetric for some reason
geog_dist <- sp::spDists(x = latlong, y = latlong, longlat = TRUE)

######################################################################################
# run conStuct analyses for Ks 1-7 (North paleo-island)
# A. d. ravitergum, A. d. ignigularis, A. d. properus, and A. d. dominiciensis lineages
# and Ks 1-5 (South paleo-island)
# A. d. dominicensis III, A. d. favillarum, A. d. aurifer, A. d. vinosus, A. d. suppar
#######################################################################################

## North paleo-island
dir.create("north-island/outputs/")

#cl <- makeCluster(10, type = "FORK")
#registerDoParallel(cl) ## Not necessary for parallelization?
myRun <- x.validation(train.prop = 0.8,
    n.reps = 8,
    K = 1:7,
    freqs = freqs,
    coords = latlong,
    geoDist = geog_dist,
    n.iter = 1e4,
    make.figs = FALSE,
    save.files = FALSE,
    parallel = TRUE,
    n.nodes = 16,
    prefix = "north-island/outputs/north-subsp-")
print(warnings())

## South paleo-island
dir.create("south-island/outputs/")

myRun <- x.validation(train.prop = 0.8,
    n.reps = 8,
    K = 1:5,
    freqs = freqs,
    coords = latlong,
    geoDist = geog_dist,
    n.iter = 1e4,
    make.figs = FALSE,
    save.files = FALSE,
    parallel = TRUE,
    n.nodes = 16,
    prefix = "south-island/outputs/south-subsp-")
print(warnings())
