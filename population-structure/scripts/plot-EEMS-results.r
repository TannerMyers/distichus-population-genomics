setwd('/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/eems/AnoDist_aln')

library(rEEMSplots)
library(rgdal)
library(rgeos)
library(raster)

border <- readOGR("/home/tcm0036/distichus-spatial-variation/data/shape-files/Hispaniola.shp")

# Create vector to loop over
demeNumbers <- c(200, 400)

for(sp in c("distichus", "brevirostris")){
    for(i in demeNumbers){
        if(sp == "distichus"){
            MCMCPATH = c(paste0(sp, '/data/distichus-subsp-EEMS-nDemes', i, '-chain1'),
                        paste0(sp, '/data/distichus-subsp-EEMS-nDemes', i, '-chain2'),
                        paste0(sp, '/data/distichus-subsp-EEMS-nDemes', i, '-chain3'))

            PLOTPATH = paste0(sp, '/plots/Demes-', i)

            rEEMSplots::eems.plots(mcmcpath = MCMCPATH, plotpath = PLOTPATH,
                    longlat = TRUE, out.png = FALSE,
                    add.grid = FALSE, add.outline = TRUE,
                    m.plot.xy = {plot(border, col = NA, add = TRUE)},
                    q.plot.xy = {plot(border, col = NA, add = TRUE)})
        } else {
            MCMCPATH = c(paste0(sp, '/data/brevirostris-EEMS-nDemes', i, '-chain1'),
                        paste0(sp, '/data//brevirostris-EEMS-nDemes', i, '-chain2'),
                        paste0(sp, '/data//brevirostris-EEMS-nDemes', i, '-chain3'))

            PLOTPATH = paste0(sp, '/plots/Demes-', i)

            rEEMSplots::eems.plots(mcmcpath = MCMCPATH, plotpath = PLOTPATH,
                    longlat = TRUE, out.png = FALSE,
                    add.grid = FALSE, add.outline = TRUE,
                    m.plot.xy = {plot(border, col = NA, add = TRUE)},
                    q.plot.xy = {plot(border, col = NA, add = TRUE)})
        }
    }
}