setwd('/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/eems/distichus-subsp/data')

library(rEEMSplots)
library(rgdal)
library(rgeos)
library(raster)
library(reemsplots2)

DOM <- getData('GADM', country='DOM', level=0)
HTI <- getData('GADM', country= 'HTI', level =0)

row.names(DOM) <- paste("DOM", row.names(DOM), sep="_")
row.names(HTI) <- paste("HTI", row.names(HTI), sep="_")

countryinfo <- rbind(HTI, DOM, makeUniqueIDs = TRUE)
border <- gSimplify(countryinfo, tol=0.01, topologyPreserve=TRUE)


demeNumbers <- c(400,600,800,1200)              # Create vector to loop over

for(i in demeNumbers){

    MCMCPATH = c(paste0('Hispaniola-distichus-ssp-pruned-EEMS-nDemes',i,'-chain1'), 
                paste0('Hispaniola-distichus-ssp-pruned-EEMS-nDemes',i,'-chain2'), 
                paste0('Hispaniola-distichus-ssp-pruned-EEMS-nDemes',i,'-chain3'))

    PLOTPATH = paste0('plots/Hispaniola-distichus-subsp-EEMS-nDemes', i,'-shapefile')


    rEEMSplots::eems.plots(mcmcpath = MCMCPATH, plotpath = PLOTPATH, 
            longlat = TRUE, out.png = FALSE, add.grid = TRUE, add.outline = TRUE,
	        m.plot.xy = { plot(border, col = NA, add = TRUE)},
	        q.plot.xy = { plot(border, col = NA, add = TRUE)})

}

 #reemsplots2::make_eems_plots(mcmcpath = MCMCPATH, longlat = TRUE)
