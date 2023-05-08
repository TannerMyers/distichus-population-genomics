setwd("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/population-structure/AnoDist/geno100")

library(tidyverse)
library(LEA)
library(sf)
library(terra)
library(tidyterra)
library(ggmap)
library(marmap)
library(scatterpie)
library(ggnewscale)

#################################################################################

## Load shapefile
#border <- rgdal::readOGR("~/Dropbox/Distichus_Project/distichus-spatial-project/data/shape-files/Hispaniola.shp")
border <- st_read(dsn = "~/Dropbox/Distichus_Project/distichus-spatial-project/data/shape-files/")

alt <- rast("~/Dropbox/Distichus_Project/distichus-spatial-project/data/elevation_new/SRTM_elevation_1km.asc")
    alt <- alt <- mask(alt, border)
alt_df <- as.data.frame(alt, xy = T)
    colnames(alt_df) <- c("x", "y", "elevation")

slope <- terrain(alt, "slope", unit = "radians")
aspect <- terrain(alt, "aspect", unit = "radians")
hill <- shade(slope, aspect, 45, 315)
    names(hill) <- "shades"
hill_df <- as.data.frame(hill, xy = TRUE)

pal_grays <- hcl.colors(1000, "Grays")
index <- hill %>%
    mutate(index_col = scales::rescale(shades, to = c(1, length(pal_grays)))) %>%
    mutate(index_col = round(index_col)) %>%
    pull(index_col)
vector_cols <- pal_grays[index]

hill_plot <- ggplot() + geom_spatraster(data = hill, fill = vector_cols, maxcell = Inf, alpha = 1, show.legend = FALSE)

basemap <- hill_plot +
    geom_raster(data = alt_df, aes(x, y, fill = elevation), alpha = 0.7, show.legend = FALSE) +
    scale_fill_hypso_tint_c(palette = "colombia", breaks = c(seq(0, 500, 50), seq(500, 2000, 500), 3000)) +
    guides(fill = guide_colorsteps(barwidth = 20, barheight = .5)) +
    theme_void() + theme(axis.line = element_blank(),
    		axis.text.x = element_blank(),
    		axis.text.y = element_blank(),
    		axis.ticks = element_blank(),
    		axis.title.x = element_blank(),
    		axis.title.y = element_blank(),
    		legend.position = "none",
   		  	panel.grid.major = element_line(color = "white", size = 0.2),
            panel.grid.minor = element_blank(),
    		plot.title = element_text(size=18, color="grey20", hjust=1, vjust=-5),
    		plot.caption = element_text(size=8, color="grey70", hjust=.15, vjust=20),
    		plot.margin = unit(c(t=0, r=0, b=0, l=0), "lines"), #added these narrower margins to enlarge maps
    		plot.background = element_rect(fill = "white", color = NA),
    		panel.background = element_rect(fill = "white", color = NA),
    		panel.border = element_blank())

## Use Hispaniola shapefile and elevation raster to generate a topographic map in gggplot
# basemap <- ggplot() + geom_tile(data = alt_df, aes(x = x, y = y, fill = elevation)) + scale_fill_etopo() +
#     coord_sf(crs = "+proj=longlat +datum=WGS84 +no_defs") +
#     geom_sf(data = border, fill = NA, color = "black") +
#     theme_minimal() + theme(text = element_text(family = "georg", color = "#22211d"),
#     		axis.line = element_blank(),
#     		axis.text.x = element_blank(),
#     		axis.text.y = element_blank(),
#     		axis.ticks = element_blank(),
#     		axis.title.x = element_blank(),
#     		axis.title.y = element_blank(),
#     		legend.position = "none",
#    		  	panel.grid.major = element_line(color = "white", size = 0.2),
#             panel.grid.minor = element_blank(),
#     		plot.title = element_text(size=18, color="grey20", hjust=1, vjust=-5),
#     		plot.caption = element_text(size=8, color="grey70", hjust=.15, vjust=20),
#     		plot.margin = unit(c(t=0, r=0, b=0, l=0), "lines"), #added these narrower margins to enlarge maps
#     		plot.background = element_rect(fill = "white", color = NA),
#     		panel.background = element_rect(fill = "white", color = NA),
#     		panel.border = element_blank()) + geom_point()

#################################################################################

# For optimal value of K, plot ancestry coefficients
    RAD_data_north <- read_table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/info/distichus-north-island-popmap-master.tsv")
    RAD_data_south <- read_table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/info/distichus-south-island-popmap-master.tsv")
    
north_snmf <- load.snmfProject("lea/distichus/north/distichus-north-island-no-miss-50kbLD.snmfProject")
    plot(north_snmf, cex = 1.2, col = "lightblue", pch = 19)
    for (i in 1:6){
        print(mean(cross.entropy(north_snmf, K = i)))
    }
    k <- 4
    ce <- cross.entropy(north_snmf, K = k)
    best_run <- which.min(ce)
        print(paste0("The lowest cross entropy score for a K of ", k, " is ", ce[best_run]))

    n_qmatrix <- LEA::Q(north_snmf, run = best_run, K = k)
    Qmatrix  <- tess3r::as.qmatrix(n_qmatrix)

    # get cluster assignments per individual
    cluster <- apply(n_qmatrix, 1, which.max)

    n_df <- as_tibble(cbind(RAD_data_north[,c("Sample_ID_pop", "Locality", "Longitude", "Latitude")], cluster, n_qmatrix))
        colnames(n_df) <- c("specimen", "locality", "longitude", "latitude", "cluster", "V1", "V2", "V3", "V4")
        # write_delim(x = n_df, file = "lea/distichus/north/north-island-subsp-K4.tsv", col_names = TRUE, delim = "\t")
        n_df <- n_df %>% arrange(cluster, longitude)

    n_colors <- c("sienna", "goldenrod1", "seashell4", "orangered")
    pdf("lea/distichus/north/north-island-subsp-K4-STRUCTURE-plot.pdf", width = 17, height = 9)
        barplot(t(n_df[, 6:ncol(n_df)]), col = n_colors, names.arg = n_df$specimen, las = 2, cex.names = 0.3, border = NA)
    dev.off()

    clusters <- grep("V", names(n_df))
    n_avg_admix <- stats::aggregate(n_df[, clusters], by = n_df[, c("locality", "longitude", "latitude")], mean)

north <- basemap + new_scale_fill() +
        geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = n_avg_admix, cols = colnames(n_avg_admix[, 4:ncol(n_avg_admix)]), alpha = 1) +
        scale_fill_manual(values = n_colors)
##############################################################################


south_snmf <- load.snmfProject("lea/distichus/south/distichus-south-island-no-miss-50kbLD.snmfProject")
    plot(south_snmf, cex = 1.2, col = "lightblue", pch = 19)
    for (i in 1:7){
        print(mean(cross.entropy(south_snmf, K = i)))
        }
    k <- 3
    ce <- cross.entropy(south_snmf, K = k)
    best_run <- which.min(ce)
        print(paste0("The lowest cross entropy score for a K of ", k, " is ", ce[best_run]))

    s_qmatrix <- LEA::Q(south_snmf, run = best_run, K = k)
    Qmatrix  <- tess3r::as.qmatrix(s_qmatrix)

    # get cluster assignments per individual
    cluster <- apply(s_qmatrix, 1, which.max)

    s_df <- as_tibble(cbind(RAD_data_south[,c("Sample_ID_pop", "Locality", "Longitude", "Latitude")], cluster, s_qmatrix))
        colnames(s_df) <- c("specimen", "locality", "longitude", "latitude", "cluster", "V1", "V2", "V3")
        write_delim(x = s_df, file = "lea/distichus/south/south-island-subsp-K3.tsv", col_names = TRUE, delim = "\t")
        s_df <- s_df %>% arrange(cluster, longitude)

    s_colors <- c("dodgerblue", "burlywood", "darkorchid2")
    pdf("lea/distichus/south/south-island-subsp-K3-STRUCTURE-plot.pdf", width = 8, height = 9)
        barplot(t(s_df[, 6:ncol(s_df)]), col = s_colors, names.arg = s_df$specimen, las = 2, cex.names = 0.3, border = NA)
    dev.off()

    clusters <- grep("V", names(s_df))
    s_avg_admix <- aggregate(s_df[, clusters], by = s_df[, c("locality", "longitude", "latitude")], mean)

north + new_scale_fill() +
    geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = s_avg_admix, cols = colnames(s_avg_admix[, 4:ncol(s_avg_admix)]), alpha = 1) +
    scale_fill_manual(values = s_colors)
ggsave("snmf-pie-plot-hierarchical.pdf")