setwd('/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/eems/distichus-subsp/data')

library(rEEMSplots)
library(rgdal)
library(rworldmap)
library(reemsplots2)

map_world <-getMap()

map_hispaniola <- map_world[c(which(map_world@data$SOVEREIGNT == "Dominican Republic"), which(map_world@data$SOVEREIGNT == "Haiti")), ]

demeNumbers <- c(400,600,800,1200)              # Create vector to loop over

for(i in demeNumbers){

    MCMCPATH = c(paste0('Hispaniola-distichus-ssp-pruned-EEMS-nDemes',i,'-chain1'), 
                paste0('Hispaniola-distichus-ssp-pruned-EEMS-nDemes',i,'-chain2'), 
                paste0('Hispaniola-distichus-ssp-pruned-EEMS-nDemes',i,'-chain3'))

    PLOTPATH = paste0('plots/Hispaniola-distichus-subsp-EEMS-nDemes', i,'-shapefile')


    rEEMSplots::eems.plots(mcmcpath = MCMCPATH, plotpath = PLOTPATH, 
            longlat = TRUE, out.png = FALSE, add.grid = TRUE, add.outline = TRUE,
	        m.plot.xy = { plot(map_hispaniola, col = NA, add = TRUE)},
	        q.plot.xy = { plot(map_hispaniola, col = NA, add = TRUE)})

}

 #reemsplots2::make_eems_plots(mcmcpath = MCMCPATH, longlat = TRUE)
