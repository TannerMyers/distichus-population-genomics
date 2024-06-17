setwd('/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/population-structure/eems/')

work_dir <- getwd()

library(rEEMSplots)
library(rworldmap)
library(rworldxtra)
library(terra)
library(sf)
library(tidyverse)

# Load points
coords <- read_table('data/distichus_geno100_filtered_LD_50kb_100_0.1_ready.coord')

projection.in = "+proj=longlat +datum=WGS84"
projection.out = "+proj=longlat +datum=WGS84"

# Create vector to loop over
demeNumbers <- c(100, 200, 400)

for(i in demeNumbers){
        MCMCPATH = c(paste0('distichus_geno100_filtered_LD_50kb_100_0.1_ready-EEMS-nDemes', i, '-chain1'),
                    paste0('distichus_geno100_filtered_LD_50kb_100_0.1_ready-EEMS-nDemes', i, '-chain2'),
                    paste0('distichus_geno100_filtered_LD_50kb_100_0.1_ready-EEMS-nDemes', i, '-chain3'))

        PLOTPATH = paste0('plots/Demes-', i)

        rEEMSplots::eems.plots(mcmcpath = MCMCPATH, plotpath = PLOTPATH,
                    longlat = TRUE, out.png = FALSE,
                    add.grid = FALSE, add.outline = TRUE,
                    m.plot.xy = {points(coords, pch = 19, cex = 0.75)},
                    q.plot.xy = {points(coords, pch = 19, cex = 0.75)},
                    add.map = TRUE, projection.in = projection.in,
                    projection.out = projection.out
                    #m.plot.xy = {plot(border, col = NA, add = TRUE)},
                    #q.plot.xy = {plot(border, col = NA, add = TRUE)}
    )
}