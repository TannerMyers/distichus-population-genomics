setwd("/mmfs1/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/gea/rda")

library(vcfR)
library(tidyverse)
library(adegenet)
library(vegan)
library(permute)
library(lattice)
library(adespatial)
library(sp)

load_funcs <- function(){
  # CA.newr
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/CA.newr.R')
  
  # PCA.newr
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/PCA.newr.R')
  
  # Rao
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/Rao.R')
  
  # bartlett.perm
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/bartlett.perm.R')
  
  # boxplerk
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/boxplerk.R')
  
  # boxplert
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/boxplert.R')
  
  # cleanplot.pca
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/cleanplot.pca.R')
  
  # coldiss
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/coldiss.R')
  
  # drawmap
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/drawmap.R')
  
  # drawmap3
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/drawmap3.R')
  
  # hcoplot
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/hcoplot.R')
  
  # panelutils
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/panelutils.R')
  
  # plot.lda
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/plot.lda.R')
  
  # plot.links
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/plot.links.R')
  
  # polyvars
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/polyvars.R')
  
  # quickMEM
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/quickMEM.R')
  
  # scalog
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/scalog.R')
  
  # screestick
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/screestick.R')
  
  # sr.value
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/sr.value.R')
  
  # triplot.rda
  source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/triplot.rda.R')
  
}
load_funcs()

## Load imputed VCF (139,964 SNPs for 191 specimens)
vcf <- read.vcfR('/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/vcf/distichus-subsp-only/imputation/distichus_geno0.8_ind0.8_maf0.05_dp4_imputed.vcf')
vcf

## Load specimen information
columns <- c('specimen', 'specimen_pop', 'filepath', 'subspecies', 'population', 'complex', 'island', 'locality', 'longitude', 'latitude', 'note')
RAD_data <- read_table("~/distichus-ddRAD/info/distichus-species-group-Hispaniola-only-popmap-cleaned.tsv", col_names = FALSE)
        colnames(RAD_data) <- columns

## Isolate per-sampling locality information
coords <- distinct(RAD_data[,c('locality', 'longitude', 'latitude')])
        dim(coords)

dat <- vcfR2genind(vcf)
        pop(dat) <- RAD_data$locality

pop_data <- genind2genpop(dat)
pop_data

print(rownames(pop_data@tab) == coords$locality) # populations ordered correctly

## Need to go from counts to frequencies
freq <- makefreq(pop_data)
        ncol(freq)
        nrow(freq)

s <- seq(from = 1, to = ncol(freq), by = 2)
freq_l1 <- freq[, s]

## Perform Hellinger transformation of allele frequencies
freq_h <- decostand(freq_l1, method = "hellinger", MARGIN = 2)

# Project coordinates
latlong <- SpatialPoints(coords[, c("longitude", "latitude")], proj4string = CRS("+proj=longlat"))

dist_xy <- as.data.frame(spTransform(latlong, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")))

# Estimate dbMEMs on a per-locality basis
dbmem_all <- dbmem(dist_xy)
        save(dbmem_all, file = 'distichus-dbMEMs.RData')

## Run quickMEM function
# dbMEM_1 <- quickMEM(freq_l1, dist_xy, perm.max = 999, method = 'fwd')
#        save(dbMEM_1, file = "non-transformed-quickMEM.RData")
dbMEM_2 <- quickMEM(freq_h, dist_xy, perm.max = 999, method = 'fwd')
        save(dbMEM_2, file = "hellinger-transformed-quickMEM.RData")