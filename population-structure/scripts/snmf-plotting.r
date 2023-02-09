setwd("Dropbox/Distichus_Project")

library(sp)
library(DescTools)
library(gridExtra)
library(raster)
library(tiff)
library(TeachingDemos)
library(rgdal)
library(RNetCDF)
library(maps)
library(mapdata)
library(shapefiles)
library(ncdf4)
library(maptools)
library(rgeos)
library(beepr)
library(ggmap)
library(GISTools)
library(ggplot2)
library(graphics)
library(prettymapr)
library(extrafont)
library(tidyverse)
library(LEA)

#################################################################################
# Make base map
bathy_fname <- "~/Dropbox/Distichus_Project/Mapping/distichus_hispaniola_map.nc" # from https://download.gebco.net/

nc <- open.nc(bathy_fname)
print.nc(nc)
#plot(nc)

tmp <- read.nc(nc)
print.nc(nc)
y <- tmp$lat
y <- as.matrix(y, nrow=length(y), ncol=1)
y<- y[order(y[,1], decreasing = F),]
y <- as.matrix(y, nrow=length(y), ncol=1)
x <- tmp$lon
x <- as.matrix(x, nrow=length(x), ncol=1)
z <- tmp$Elevation.relative.to.sea.level #layer
rm(tmp)
close.nc(nc)

xlim <- range(x)
ylim <- range(y)


###Plot
#make palette
#ocean.pal <- colorRampPalette(c("#000000", "#000209", "#000413", "#00061E", "#000728", "#000932", "#002650", "#00426E", "#005E8C", "#007AAA", "#0096C8", "#22A9C2", "#45BCBB", "#67CFB5", "#8AE2AE", "#ACF6A8", "#BCF8B9", "#CBF9CA", "#DBFBDC", "#EBFDED")) 

ocean.pal <- colorRampPalette(c("#064273","#1da2d8", "#def3f6", "#76b6c4","#7fcdff")) 

land.pal <- colorRampPalette(
  c("#336600", "#F3CA89", "#D9A627", 
    "#A49019", "#9F7B0D", "#996600", "#B27676", "#C2B0B0", "#E5E5E5", 
    "#FFFFFF")
)

#New land colors
#land.pal <- colorRampPalette(c("#6F9C57","#A2BB74","#D9DF92","#F2E8A2","#F8DB93","#F3CE86","#DEAE84", "#C59781","#D29C97", "#E9C3D0","#FEF3F7"))

#land.pal <- colorRampPalette(c("grey34","grey81"))
zbreaks <- seq(-20000, 6000, by=10)
cols <-c(ocean.pal(sum(zbreaks<=0)-1), land.pal(sum(zbreaks>0)))

## Load shapefile
DOM <- getData('GADM', country = 'DOM', level=0, path = paste0("~/Dropbox/Distichus_Project/distichus-spatial-project/data/shape-files/"), download = F)
HTI <- getData('GADM', country = 'HTI', level=0, path = paste0("~/Dropbox/Distichus_Project/distichus-spatial-project/data/shape-files/"), download = F)

row.names(DOM) <- paste("DOM", row.names(DOM), sep = "_")
row.names(HTI) <- paste("HTI", row.names(HTI), sep = "_")

countryinfo <- rbind(HTI, DOM, makeUniqueIDs = TRUE)
border <- gSimplify(countryinfo, tol=0.01, topologyPreserve=TRUE)
#################################################################################

# For optimal value of K, plot ancestry coefficients
for (i in 5:10){
    k = i
    ce <- cross.entropy(obj.snmf, K = i) # K = 10 determined optimal
    best_run <- which.min(ce)
        print(paste0("The lowest cross entropy score for a K of ", i, " is ", ce[best_run]))

    qmatrix <- LEA::Q(obj.snmf, run = best_run, K = k)
    Qmatrix  <- tess3r::as.qmatrix(qmatrix)

    # get cluster assignments per individual
    cluster <- apply(qmatrix, 1, which.max)

    # Merge Qmatrix with specimen information
    sNMF_ancestry_coefficients_cluster <- as.data.frame(cbind(RAD_data3[,c("Sample_ID_pop", "Locality", "Longitude", "Latitude")], cluster, qmatrix))
    ## Order by cluster
    snmf <- sNMF_ancestry_coefficients_cluster[order(as.numeric(sNMF_ancestry_coefficients_cluster$cluster)), ]
        write_delim(x = snmf,
            file = paste0("sNMF_K", k, "_ancestry_coefficients_cluster-specimen-cleaned.tsv"),
            col_names = TRUE, delim = "\t")

    clusters <- grep("V", names(snmf))
    
    avg_admix <- aggregate(snmf[, clusters], by = snmf[, c("Locality", "Longitude", "Latitude")], mean)
    #    avg_admix <- avg_admix[order(as.character(avg_admix$Locality)), ]
    
    png(paste0("snmf-pie-plot-K", k, ".pdf"))
        raster::image(x,y,z=as.matrix(z), col=cols, breaks=zbreaks, useRaster = T,xlab = "", ylab = "", axes = F)
        plot(border, add=T, lwd=0.2, border="black")
        box(lwd=1.5)
        scalebar(type = "bar", divs = 2, d = 200, below = "km", xy = NULL)
        addnortharrow(scale = 0.35, pos = "topright")
            for (i in 1:length(unique(avg_admix$Locality))){
                mapplots::add.pie(z = as.numeric(avg_admix[i, 4:ncol(avg_admix)]),
                                x = avg_admix[i, 2], y = avg_admix[i, 3],
                                labels = "", col = paste0("cols_", k), radius = 0.07)
            }
        dev.off()
}

snmf.k.5 <- read_table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/population-structure/lea/sNMF/distichus-subsp/sNMF_K5_ancestry_coefficients_cluster-specimen-cleaned.tsv", 
                    col_names = T)
clusters <- grep("V", names(snmf.k.5))
avg_admix <- aggregate(snmf.k.5[, clusters], by = snmf.k.5[, c("Locality", "Longitude", "Latitude")], mean)
    avg_admix <- avg_admix[order(as.character(avg_admix$Locality)), ]

cols_5 <- c("#634102", "#1b0e6e", "#cf5d34", "#fddc1f", "gray80")

pdf("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/population-structure/lea/sNMF/distichus-subsp/snmf-K5-pie-map.pdf",
    width = 8, height = 5)
    par(mar=c(1,1,1,1))
    raster::image(x,y,z=as.matrix(z), col=cols, breaks=zbreaks, useRaster=T,xlab="", ylab="",axes=F)
    plot(border, add=T, lwd=0.2, border="black")
    box(lwd=1.5)
    scalebar(type = "bar", divs = 2, d = 200, below = "km", xy=NULL)
    addnortharrow(scale = 0.35, pos = "topright")
        for (i in 1:length(unique(avg_admix$Locality))){
            mapplots::add.pie(z = as.numeric(avg_admix[i, 4:ncol(avg_admix)]),
                            x = avg_admix[i, 2], y = avg_admix[i, 3],
                            labels = "", col = cols_5, radius = 0.05)
        }
dev.off()