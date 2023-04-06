setwd("/mmfs1/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/conStruct/AnoDist_aln/distichus-subsp")

## Here's the code used to generate files in STRUCTURE format for easy
## conversion to conStruct inputs:

# plink --allow-extra-chr --double-id --out north-island-subsp --recode structure --set-missing-var-ids @:# --vcf /scratch/tcm0036/distichus-ddRAD/analyses/population-structure/eems/AnoDist_aln/distichus/distichus-north-island-no-miss-50kbLD.recode.vcf
# plink --allow-extra-chr --double-id --out south-island-subsp --recode structure --set-missing-var-ids @:# --vcf /scratch/tcm0036/distichus-ddRAD/analyses/population-structure/eems/AnoDist_aln/distichus/distichus-south-island-no-miss-50kbLD.recode.vcf


library(tidyverse)
library(geosphere)
library(raster)
library(conStruct)
library(doParallel)
library(foreach)

print(Sys.time())

#########
## Generate input files for conStruct
#########

## (1) Load allele frequencies file
## First, make allele freq file for whole island
    # structure2conStruct(infile = "distichus-subsp.recode.strct_in",
    #     onerowperind = TRUE,
    #     start.loci = 3,
    #     start.samples = 3,
    #     missing.datum = 0, # This data matrix lacks any missing data so I just chose 0
    #     outfile = "distichus-conStruct-input")
# load("distichus-conStruct-input.RData")

## Now, for north paleo-island subset
structure2conStruct(infile = "north-island/north-island-subsp.recode.strct_in",
	 onerowperind = TRUE,
	 start.loci = 3,
	 start.samples = 3,
 	 missing.datum = 0,
	 outfile = "north-island/north-subsp-input")
load("north-island/north-subsp-input.RData")

## Last, for south paleo-island subset
structure2conStruct(infile = "south-island/south-island-subsp.recode.strct_in",
	 onerowperind = TRUE,
	 start.loci = 3,
	 start.samples = 3,
	 missing.datum = 0,
	 outfile = "south-island/south-subsp-input")
load("south-island/south-subsp-input.RData")

# (2) Load coordinates for specimen sampling localities

RAD_data <- read_table("~/distichus-ddRAD/info/distichus-subsp-Hispaniola-popmap-cleaned.tsv", col_names = TRUE)
longlat <- as.matrix(RAD_data[,c('Longitude', 'Latitude')])

## North paleo-island pop map
RAD_data <- read_table("~/distichus-ddRAD/info/distichus-north-island-popmap-master.tsv", col_names = TRUE)
longlat <- as.matrix(RAD_data[,c('Longitude', 'Latitude')])

## South paleo-island pop map
RAD_data <- read_table("~/distichus-ddRAD/info/distichus-south-island-popmap-master.tsv", col_names = TRUE)
longlat <- as.matrix(RAD_data[,c('Longitude', 'Latitude')])

## (3) Generate geographic distance matrix from specimen sampling localities

#geog_dist <- sp::spDists(x = longlat, y = longlat, longlat = TRUE)
geog_dist <- geosphere::distm(x = longlat, y = longlat, fun = distGeo)

##############################################################################################################
# run conStuct analyses for Ks 1-5 (North paleo-island)
# A. d. ravitergum, A. d. ignigularis, A. d. properus, and A. d. dominiciensis II lineage
# and Ks 1-4 (South paleo-island)
# A. d. dominicensis I & IV, A. d. favillarum, and A. d. dom III, A. d. aurifer, A. d. vinosus, A. d. suppar
##############################################################################################################

# sbatch --job-name=north-conStruct --partition=jro0014_amd --mem=100G --time=27-24:00:00 --ntasks=16 --wrap="Rscript run-conStruct-distichus.r"
# sbatch --job-name=south-conStruct --partition=jro0014_amd --mem=100G --time=14-24:00:00 --ntasks=16 --wrap="Rscript run-conStruct-distichus.r"

## North paleo-island
# dir.create("north-island/outputs/")

#cl <- makeCluster(10, type = "FORK")
#registerDoParallel(cl) ## Not necessary for parallelization?
myRun <- x.validation(train.prop = 0.8,
    n.reps = 8,
    K = 1:5,
    freqs = freqs,
    coords = longlat,
    geoDist = geog_dist,
    n.iter = 1e4,
    make.figs = TRUE,
    save.files = TRUE,
    parallel = TRUE,
    n.nodes = 16,
    prefix = "north-island/outputs/north-subsp")
print(warnings())

## South paleo-island
# dir.create("south-island/outputs/")

myRun <- x.validation(train.prop = 0.8,
    n.reps = 8,
    K = 1:4,
    freqs = freqs,
    coords = longlat,
    geoDist = geog_dist,
    n.iter = 1e4,
    make.figs = TRUE,
    save.files = TRUE,
    parallel = TRUE,
    n.nodes = 16,
    prefix = "south-island/outputs/south-subsp")
print(warnings())