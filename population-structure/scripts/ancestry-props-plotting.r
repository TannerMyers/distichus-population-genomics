## This script infer population structure using both multivariate and model-based approaches
## and plots ancestry proportions inferred by population structure approaches
## snmf and ADMIXTURE as both barplots and pie-charts overlaid onto a map of Hispaniola
###############################################################################################

## Set working directory
work_dir <- "~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/AnoDist_aln/population-structure/ADMIXTURE_PCA" #CHANGE IF NEEDED; getwd()
setwd(work_dir)

## Load necessary packages
library(tidyverse)
library(vcfR)
library(LEA)
library(adegenet)
library(sf)
library(terra)
library(tidyterra)
library(marmap)
library(scatterpie)
library(ggnewscale)
library(cowplot)

################################################################
## Load specimen information into R
columns <- c("specimen", "specimen_pop", "filepath", "subspecies", "population", "complex", "island", "locality", "longitude", "latitude", "note")
RAD_data <- read_table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/info/distichus-species-group-Hispaniola-only-popmap-cleaned.tsv", col_names = FALSE) 
RAD_data_north <- read_table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/info/north-island-dataset-master.tsv", col_names = FALSE)
RAD_data_south <- read_table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/info/south-island-dataset-master.tsv", col_names = FALSE)
    colnames(RAD_data) <- columns
    colnames(RAD_data_north) <- columns
    colnames(RAD_data_south) <- columns

## Load SNP datasets into R: MISSINGNESS
vcf_100 <- read.vcfR("../../distichus_geno100_filtered_LD_50kb_100_0.1_biallelic_pruned.recode.vcf")
    colnames(vcf_100@gt)[2:192] <- RAD_data$specimen_pop
vcf_90 <- read.vcfR("../../distichus_geno90_filtered_LD_50kb_100_0.1_biallelic_pruned.recode.vcf")
    colnames(vcf_90@gt)[2:192] <- RAD_data$specimen_pop
vcf_80 <- read.vcfR("../../distichus_geno80_filtered_LD_50kb_100_0.1_biallelic_pruned.recode.vcf")
    colnames(vcf_80@gt)[2:192] <- RAD_data$specimen_pop

## Load SNP datasets into R: HIERARCHICAL GEOGRAPHIC PARTITIONS
#n_vcf <- read.vcfR("../north_geno100_50kb_100_0.1_biallelic.recode.vcf")
n_vcf <- read.vcfR("../north_geno90_50kb_100_0.1_biallelic.recode.vcf")
    colnames(n_vcf@gt)[2:144] <- RAD_data_north$specimen_pop
#s_vcf <- read.vcfR("../south_geno100_50_100_0.1_biallelic.recode.vcf")
s_vcf <- read.vcfR("../south_geno90_50kb_100_0.1_biallelic.recode.vcf")
    colnames(s_vcf@gt)[2:49] <- RAD_data_south$specimen_pop

### Convert to adegenet objects
## 100%
gen_100 <- vcfR2genlight(vcf_100)
    pop(gen_100) <- RAD_data$subspecies
    indNames(gen_100) <- RAD_data$specimen_pop
## 90%
gen_90 <- vcfR2genlight(vcf_90)
    pop(gen_90) <- RAD_data$subspecies
## 80
gen_80 <- vcfR2genlight(vcf_80)
    pop(gen_80) <- RAD_data$subspecies

## north island (100%)
n_gen <- vcfR2genlight(n_vcf)
    pop(n_gen) <- RAD_data_north$subspecies

## south island (100%)
s_gen <- vcfR2genlight(s_vcf)
    pop(s_gen) <- RAD_data_south$subspecies

#################################################################################
## (1) Visualize population structure with multivariate approach (PCA)
#################################################################################

# This performs PCA on the genlight object `gen` and requires manually choosing the number of principal components to retain
## I selected 4 axes to retain for all 3 datasets
pca_100 <- glPca(gen_100, center = TRUE, scale = FALSE, nf = 4)
pca_90 <- glPca(gen_90, center = TRUE, scale = FALSE, nf = 4)
pca_80 <- glPca(gen_80, center = TRUE, scale = FALSE, nf = 4)

pc100_df <- as_tibble(cbind(gen_100@ind.names, RAD_data$subspecies, as.data.frame(pca_100$scores)))
pc90_df <- as_tibble(cbind(gen_90@ind.names, RAD_data$subspecies, as.data.frame(pca_90$scores)))
pc80_df <- as_tibble(cbind(gen_80@ind.names, RAD_data$subspecies, as.data.frame(pca_80$scores)))
    colnames(pc100_df) <- c("specimen", "population", "PC1", "PC2", "PC3", "PC4")
    colnames(pc90_df) <- c("specimen", "population", "PC1", "PC2", "PC3", "PC4")
    colnames(pc80_df) <- c("specimen", "population", "PC1", "PC2", "PC3", "PC4")

barplot(100*pca_100$eig/sum(pca_100$eig), col = heat.colors(100), main="PCA (100%) Eigenvalues")
barplot(100*pca_90$eig/sum(pca_90$eig), col = heat.colors(50), main="PCA (90%) Eigenvalues")
barplot(100*pca_80$eig/sum(pca_80$eig), col = heat.colors(50), main="PCA (80%) Eigenvalues")

## Set color palette
myColors <- c("dodgerblue1", "burlywood", "sienna", "dodgerblue2", "burlywood1", "darkorchid2", "orangered", "gray80", "yellow", "dodgerblue3", "dodgerblue4")
col_labs <- c("A. d. aurifer", "A. d. dominicensis I", "A. d. dominicensis II", "A. d. dominicensis III", "A. d. dominicensis IV",
            "A. d. favillarum", "A. d. ignigularis", "A. d. properus", "A. d. ravitergum", "A. d. suppar", "A. d. vinosus")

## 100%
plt_100 <- ggplot(data = pc100_df, aes(x = PC1, y = PC2, color = population, label = specimen, alpha = 0.85)) +
    geom_point(cex = 3) + scale_color_manual(values = myColors, labels = col_labs,
        guide = guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
    guides(alpha = "none") +
    labs(color = "Subspecies", x = "Axis 1 (14.54% of variation)", y = "Axis 2 (13.37% of variation)") + # UDPATE values if needed
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 3),
    axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"))
plt_100
ggsave("PCA-100-fulldata.pdf")

## Now, 90 %
plt_90 <- ggplot(data = pc90_df, aes(x = PC1, y = PC2, color = population, label = specimen)) +
    geom_point(cex = 5, alpha = 1, show.legend = FALSE) + # change to TRUE for legend
    scale_color_manual(values = myColors, labels = col_labs, guide = guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
    labs(x = "Axis 1 (17.33% of variation)", y = "Axis 2 (12.66% of variation)") + # UDPATE values if needed
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 3), panel.background = element_rect(fill = "white"),
        axis.title = element_text(size = 22, colour = "black", family = "sans", face = "plain"),
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16))
plt_90
plt_90 + geom_mark_ellipse(data = pc90_df, aes(color = cluster)) + new_scale_color() + scale_color_manual(values = K_5_cols)
ggsave("PCA-90-fulldata.pdf")

## Now, 80%
plt_80 <- ggplot(data = pc80_df, aes(x = PC1, y = PC2, color = population, label = specimen, alpha = 0.85)) +
    geom_point(cex = 3) + scale_color_manual(values = myColors, labels = col_labs,
        guide = guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
    guides(alpha = "none") +
    labs(color = "Subspecies", x = "Axis 1 (17.13% of variation)", y = "Axis 2 (11.54% of variation)") + # UDPATE values if needed
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 3),
    axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"))
plt_80
ggsave("PCA-80-fulldata.pdf")

## Now, do hierarchical subsets
n_pca <- glPca(n_gen, center = TRUE, scale = FALSE, nf = 3) # Argument for 3, went with 4
    barplot(100*n_pca$eig/sum(n_pca$eig), col = heat.colors(80), main="PCA Eigenvalues")
s_pca <- glPca(s_gen, center = TRUE, scale = FALSE, nf = 2) # probably 2, went with 3
    barplot(100*s_pca$eig/sum(s_pca$eig), col = heat.colors(50), main="PCA Eigenvalues")

npc_df <- as_tibble(cbind(n_gen@ind.names, RAD_data_north$subspecies, as.data.frame(n_pca$scores)))
    colnames(npc_df) <- c("specimen", "population", "PC1", "PC2", "PC3") 
spc_df <- as_tibble(cbind(s_gen@ind.names, RAD_data_south$subspecies, as.data.frame(s_pca$scores)))
    colnames(spc_df) <- c("specimen", "population", "PC1", "PC2", "PC3")

n_cols <- c("sienna", "orangered", "gray80", "yellow")
n_labs <- c("A. d. dominicensis II", "A. d. ignigularis", "A. d. properus", "A. d. ravitergum")
n_plt <- ggplot(data = npc_df, aes(x = PC1, y = PC2, color = population, label = specimen, alpha = 0.85)) +
    geom_point(cex = 3) + scale_color_manual(values = n_cols, labels = n_labs,
        guide = guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
    guides(alpha = "none") +
    labs(color = "Subspecies", x = "Axis 1 (18.58% of variation)", y = "Axis 2 (14.72% of variation)") + # UDPATE values if needed
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 3),
    axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"))
n_plt
ggsave("north-island-90complete-pc1-pc2.pdf")

s_cols <- c("dodgerblue4", "burlywood", "mediumblue", "burlywood3", "darkorchid2", "royalblue3", "royalblue")
s_labs <- c("A. d. aurifer", "A. d. dominicensis I", "A. d. dominicensis III", "A. d. dominicensis IV", "A. d. favillarum", "A. d. suppar", "A. d. vinosus")
s_plt <- ggplot(data = spc_df, aes(x = PC1, y = PC2, color = population, label = specimen, alpha = 0.85)) +
    geom_point(cex = 3) + scale_color_manual(values = s_cols, labels = s_labs,
        guide = guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
    guides(alpha = "none") +
    labs(color = "Subspecies", x = "Axis 1 (28.28% of variation)", y = "Axis 2 (17.37% of variation)") + # UDPATE values if needed
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 3),
    axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"))
s_plt
ggsave("south-island-90complete-pc1-pc2.pdf")

#################################################################################
## (2) Make hillshade map to serve as background for ancestry prop. pie-maps
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
basemap

#################################################################################
## (3) Load ADMIXTURE results into R and examine ancestry proportions
#################################################################################

## Plot CV error support
##change directory to one with ADMIXTURE log files
for (i in c(80, 90, 100)){
    #system("grep -h CV log*out > cross_validation.txt")
    cv <- read_table(paste0("geno", i, "/cross_validation.txt"), col_names = FALSE)

    ## Clean table and re-order numerically
    cv$X3 <- gsub("[\\(\\)]", "", regmatches(cv$X3, gregexpr("\\(.*?\\)", cv$X3))) %>% str_remove("K=")
        cv$X3 <- as.numeric(cv$X3)
        cv <- cv[order(cv$X3), ]

    cv_plt <- ggplot(cv, aes(x = X3, y = X4)) + geom_line() + scale_x_continuous(breaks=c(2, 4, 6, 8, 10)) +
            labs(title = "ADMIXTURE Cross-Validation Error for K=1 to K=10",
                x = "K", y = "Cross-Validation Error") +
            theme(axis.text.x = element_text(colour = "black")) +
            theme(legend.title = element_blank()) + theme(axis.text.y = element_text(colour = "black", size=12)) +
            theme(axis.text.x = element_text(colour = "black", size = 12)) +
            theme(panel.border = element_rect(colour = "black", fill = NA, size = 3),
                axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"))
    #ggsave(paste0("ADMIXTURE_CV_", i, ".pdf"))
    assign(paste0("cv_", i, "_plt"), cv_plt)
}

## Full dataset
## ADMIXTURE CV levels out around K = 4 with slight increases for K = 7 and K = 9, let"s examine how it clusters specimens
## under different numbers of ancestral populations
for (m in c(80, 90, 100)){
    for (K in c(5, 7, 8)){
        adm_Q <- read_table(paste0("geno", m, "/distichus_geno", m, "_filtered_LD_50kb_100_0.1_ready.", K, ".Q"), col_names = FALSE)
            cluster <- apply(adm_Q, 1, which.max)
            adm_Q <- as_tibble(cbind(RAD_data[, c("specimen_pop", "locality", "longitude", "latitude")], cluster, adm_Q))
            # ## Preserve order of columns by matching to order in K = 4 dataset
            # if (K == 5){
            #     adm_Q <- adm_Q %>% arrange(cluster, longitude)
            #     order <- adm_Q$specimen_pop
            # } else {
            #     adm_Q <- adm_Q[match(order, adm_Q$specimen_pop), ]
            # }
            
    bar_df <- adm_Q %>% pivot_longer(cols = starts_with("X"), names_to = "pop", values_to = "prob", names_prefix = "X")
    
    ggplot(bar_df, aes(x = reorder(specimen_pop, cluster), y = prob, fill = pop)) +
        geom_col(color = "gray", size = 0.1) +
        scale_fill_manual(values = funky(K)) +
        theme_minimal() + labs(x = "Individuals", title = paste0("K = ", K), y = "Ancestry") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_discrete(expand = expand_scale(add = 1)) +
        theme(panel.spacing.x = unit(0.1, "lines"), legend.position = "none",
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("distichus-subsp-filtered-", m, "-", K, "-barplot.pdf"))

        clusters <- grep("X", names(adm_Q))
        avg_admix <- aggregate(adm_Q[, clusters], by = adm_Q[, c("locality", "longitude", "latitude")], mean)
        ## Map pie-plots
        basemap + new_scale_fill() +
            geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = avg_admix, cols = colnames(avg_admix[, 4:ncol(avg_admix)]), alpha = 1, pie_scale = 0.7) +
            scale_fill_manual(values = funky(K))
        ggsave(paste0("distichus-subsp-filtered-", m, "-pieplot-K-", K, ".pdf"))
        
        # Keep outputs in environment
        assign(paste0("K_", K, "_geno_", m), adm_Q)
        assign(paste0("K_", K, "_geno_", m, "bar"), bar_df)
        assign(paste0("K_", K, "_geno_", m, "_admix"), avg_admix)
    }
}

## Below is for figure
K_5_cols <- c("dodgerblue4", "sienna", "yellow", "gray80", "orangered")
pie_90 <- basemap + new_scale_fill() +
            geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = K_5_geno_90_admix, cols = colnames(K_5_geno_90_admix[, 4:ncol(K_5_geno_90_admix)]), alpha = 1, pie_scale = 0.7) +
            scale_fill_manual(values = K_5_cols)
ggsave("../../../ms-figures/K_5_90_full_dataset_pieplot.pdf")
bar_90 <- ggplot(K_5_geno_90bar, aes(x = reorder(specimen_pop, cluster), y = prob, fill = pop)) +
        geom_col(color = "gray", size = 0.1) +
        scale_fill_manual(values = K_5_cols) +
        theme_minimal() + labs(x = "Individuals", title = "K = 5", y = "Ancestry") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_discrete(expand = expand_scale(add = 1)) +
        theme(panel.spacing.x = unit(0.1, "lines"), legend.position = "none",
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"),
            plot.title = element_text(hjust = 0.5))
ggsave("../../../ms-figures/K_5_90_full_dataset_barplot.pdf")
# row2 <- plot_grid(pie_90, cv_90_plt, labels = "", rel_widths = c(5,1))
# plot_grid(plt_90, row2, bar_90, labels = c("A", "B", "C"), ncol = 1, rel_heights = c(3, 2, 1))
# ggsave("test.pdf")


##################### Hierarchical analyses
# Plot CV errors
for (i in c("north", "south")){
    #system("grep -h CV log*out > cross_validation.txt")
    cv <- read_table(paste0("hierarchical-analyses/", i, "/geno90/cross_validation.txt"), col_names = FALSE)

    ## Clean table and re-order numerically
    cv$X3 <- gsub("[\\(\\)]", "", regmatches(cv$X3, gregexpr("\\(.*?\\)", cv$X3))) %>% str_remove("K=")
        cv$X3 <- as.numeric(cv$X3)
        cv <- cv[order(cv$X3), ]

    cv_plt <- ggplot(cv[1:5,], aes(x = X3, y = X4)) + geom_line() + scale_x_continuous(breaks=c(1, 3, 5)) +
            labs(title = paste0("ADMIXTURE Cross-Validation Error for K=1 to K=5"), x = "K", y = "Cross-Validation Error") +
            theme(axis.text.x = element_text(colour = "black")) +
            theme(legend.title = element_blank()) + theme(axis.text.y = element_text(colour = "black", size=12)) +
            theme(axis.text.x = element_text(colour = "black", size = 12)) +
            theme(panel.border = element_rect(colour = "black", fill = NA, size = 3),
                axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"))
    ggsave(paste0("ADMIXTURE_CV_", i, ".pdf"))
    assign(paste0(i, "_cv"), cv_plt)
}

## North paleo-island

for (K in c(3, 4, 5)){
    adm_Q <- read_table(paste0("hierarchical-analyses/north/geno90/north_geno90_50kb_100_0.1_ready.", K, ".Q"), col_names = FALSE)
        cluster <- apply(adm_Q, 1, which.max)
        adm_Q <- as_tibble(cbind(RAD_data_north[, c("specimen_pop", "locality", "longitude", "latitude")], cluster, adm_Q))
        adm_Q <- adm_Q %>% arrange(cluster)

   north_bar <- adm_Q %>% pivot_longer(cols = starts_with("X"), names_to = "pop", values_to = "prob", names_prefix = "X")
    
    ggplot(north_bar, aes(x = reorder(specimen_pop, cluster), y = prob, fill = pop)) +
        geom_col(color = "gray", size = 0.1) +
        scale_fill_manual(values = funky(K)) +
        theme_minimal() + labs(x = "Individuals", title = paste0("K = ", K), y = "Ancestry") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_discrete(expand = expand_scale(add = 1)) +
        theme(panel.spacing.x = unit(0.1, "lines"), legend.position = "none",
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("north-island-geno90-ADMIXTURE-", K, "-barplot.pdf"))

    clusters <- grep("X", names(adm_Q))
    avg_admix <- aggregate(adm_Q[, clusters], by = adm_Q[, c("locality", "longitude", "latitude")], mean)

    ## Map pie-plots
    basemap + new_scale_fill() +
        geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = avg_admix, cols = colnames(avg_admix[, 4:ncol(avg_admix)]), alpha = 1) +
        scale_fill_manual(values = funky(K))
    ggsave(paste0("north-island-ADMIXTURE-pieplot-K-", K, ".pdf"))

    # Keep outputs in environment
    assign(paste0("north_K_", K), adm_Q)
    assign(paste0("north_K_", K, "_bar"), north_bar)
    assign(paste0("north_K_", K, "_admix"), avg_admix)
}

## South paleo-island
for (K in c(3)){
    adm_Q <- read_table(paste0("hierarchical-analyses/south/geno90/south_geno90_50_100_0.1_ready.", K, ".Q"), col_names = FALSE)
        cluster <- apply(adm_Q, 1, which.max)
        adm_Q <- as_tibble(cbind(RAD_data_south[, c("specimen_pop", "locality", "longitude", "latitude")], cluster, adm_Q))
        adm_Q <- adm_Q %>% arrange(cluster)
        
       south_bar <- adm_Q %>% pivot_longer(cols = starts_with("X"), names_to = "pop", values_to = "prob", names_prefix = "X")
    
    ggplot(south_bar, aes(x = reorder(specimen_pop, cluster), y = prob, fill = pop)) +
        geom_col(color = "gray", size = 0.1) +
        scale_fill_manual(values = funky(K)) +
        theme_minimal() + labs(x = "Individuals", title = paste0("K = ", K), y = "Ancestry") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_discrete(expand = expand_scale(add = 1)) +
        theme(panel.spacing.x = unit(0.1, "lines"), legend.position = "none",
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("south-island-geno90-ADMIXTURE-", K, "-barplot.pdf"))

    clusters <- grep("X", names(adm_Q))
    avg_admix <- aggregate(adm_Q[, clusters], by = adm_Q[, c("locality", "longitude", "latitude")], mean)

    ## Map pie-plots
    basemap + new_scale_fill() +
        geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = avg_admix, cols = colnames(avg_admix[, 4:ncol(avg_admix)]), alpha = 1) +
        scale_fill_manual(values = funky(K))
    ggsave(paste0("south-island-ADMIXTURE-pieplot-K-", K, ".pdf"))

    # Keep outputs in environment
    assign(paste0("south_K_", K), adm_Q)
    assign(paste0("south_K_", K, "_bar"), south_bar)
    assign(paste0("south_K_", K, "_admix"), avg_admix)
}

######### For figure below
n_cols <- c("yellow", "sienna", "orangered", "gray80")
north_pie <- basemap + new_scale_fill() +
            geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = north_K_4_admix, cols = colnames(north_K_4_admix[, 4:ncol(north_K_4_admix)]), alpha = 1, pie_scale = 1) +
            scale_fill_manual(values = n_cols)
ggsave("../../../ms-figures/north_K_4_90_full_dataset_pieplot.pdf")
nbar_plt <- ggplot(north_K_4_bar, aes(x = reorder(specimen_pop, cluster), y = prob, fill = pop)) +
        geom_col(color = "gray", size = 0.1) +
        scale_fill_manual(values = n_cols) +
        theme_minimal() + labs(x = "Individuals", title = "K = 4", y = "Ancestry") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_discrete(expand = expand_scale(add = 1)) +
        theme(panel.spacing.x = unit(0.1, "lines"), legend.position = "none",
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"),
            plot.title = element_text(hjust = 0.5))
ggsave("../../../ms-figures/north_K_4_90_barplot.pdf")

s_cols <- c("dodgerblue", "burlywood", "darkorchid2")
south_pie <- basemap + new_scale_fill() +
            geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = south_K_3_admix, cols = colnames(south_K_3_admix[, 4:ncol(south_K_3_admix)]), alpha = 1, pie_scale = 1) +
            scale_fill_manual(values = s_cols)
ggsave("../../../ms-figures/south_K_3_90_full_dataset_pieplot.pdf")
sbar_plt <- ggplot(south_K_3_bar, aes(x = reorder(specimen_pop, cluster), y = prob, fill = pop)) +
        geom_col(color = "gray", size = 0.1) +
        scale_fill_manual(values = s_cols) +
        theme_minimal() + labs(x = "Individuals", title = "K = 3", y = "Ancestry") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_discrete(expand = expand_scale(add = 1)) +
        theme(panel.spacing.x = unit(0.1, "lines"), legend.position = "none",
            #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_text(size = 18, colour = "black", family = "Helvetica", face = "bold"),
            plot.title = element_text(hjust = 0.5))
ggsave("../../../ms-figures/south_K_3_90_barplot.pdf")


#################################################################################
## (4) Use `snmf` implemented in LEA to infer ancestry proportions
#################################################################################
## For optimal value of K, plot ancestry coefficients
## Full dataset with alpha (regularization parameter) set to 10
# snmf <- snmf()
snmf <- load.snmfProject("lea/distichus/alpha_10/distichus-subsp_variants_no-miss-50KbLD.snmfProject")
    plot(snmf, cex = 1.2, col = "lightblue", pch = 19)
    for (i in 1:10){
        print(paste("For K set to", i, "mean cross entropy is", mean(cross.entropy(snmf, K = i))))
    }

k <- 5
ce <- cross.entropy(snmf, K = k)
    best_run <- which.min(ce)
    print(paste0("The lowest cross entropy score for a K of ", k, " is ", ce[best_run]))
qmatrix <- LEA::Q(snmf, run = best_run, K = k)
cluster <- apply(qmatrix, 1, which.max)

df <- as_tibble(cbind(RAD_data[,c("Sample_ID_pop", "Locality", "Longitude", "Latitude")], cluster, qmatrix))
    colnames(df) <- c("specimen", "locality", "longitude", "latitude", "cluster", "V1", "V2", "V3", "V4", "V5")
    write_delim(x = df, file = "lea/distichus/full-dataset-K-5_alpha-10.tsv", col_names = TRUE, delim = "\t")
    clust_ord <- c(2, 4, 5, 1, 3)
    df <- df %>% arrange(cluster, longitude)
    df <- df[order(factor(df$cluster, levels = clust_ord)),]

colors <- c("orangered", "dodgerblue3", "seashell4", "goldenrod1", "sienna")
    pdf("lea/distichus/K_5_alpha_10_full_barplot.pdf", width = 17, height = 9)
        barplot(t(df[, 6:ncol(df)]), col = colors, names.arg = df$specimen, las = 2, cex.names = 0.3, border = NA)
    dev.off()

    clusters <- grep("V", names(df))
    avg_admix <- stats::aggregate(df[, clusters], by = df[, c("locality", "longitude", "latitude")], mean)

full <- basemap + new_scale_fill() +
        geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = avg_admix, cols = colnames(avg_admix[, 4:ncol(avg_admix)]), alpha = 1, pie_scale = 0.7) +
        scale_fill_manual(values = colors)
full
ggsave("K_5_alpha_10_full_dataset_non_hierarchical.pdf")

#########
## North paleo-island sans A. d. dominicensis I & IV: `east` dataset
# north_snmf <- snmf("distichus-north-island-no-miss-50kbLD.geno", K = 1:6, project = "new", repetitions = 50, tolerance = 0.00001, seed = 16, alpha = 10, entropy = TRUE, ploidy = 2)
north_snmf <- load.snmfProject("lea/distichus/north/distichus-north-island-no-miss-50kbLD.snmfProject")
    pdf("lea/distichus/north/snmf_alpha_10_K_1-6.pdf")
        plot(north_snmf, cex = 1.2, col = "lightblue", pch = 19) # Plots minimum cross entropy per value of K
    dev.off()
    for (i in 1:6){
        print(paste("For K set to", i, "mean cross entropy is", mean(cross.entropy(north_snmf, K = i))))
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
        write_delim(x = n_df, file = "lea/distichus/north/north-island-subsp-K4.tsv", col_names = TRUE, delim = "\t")
        n_df <- n_df %>% arrange(cluster, longitude)

    n_colors <- c("sienna", "goldenrod1", "seashell4", "orangered")
    pdf("lea/distichus/north/K4-alpha-10-north-island-subsp-barplot.pdf", width = 17, height = 9)
        barplot(t(n_df[, 6:ncol(n_df)]), col = n_colors, names.arg = n_df$specimen, las = 2, cex.names = 0.3, border = NA)
    dev.off()

    clusters <- grep("V", names(n_df))
    n_avg_admix <- stats::aggregate(n_df[, clusters], by = n_df[, c("locality", "longitude", "latitude")], mean)

north <- basemap + new_scale_fill() +
        geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = n_avg_admix, cols = colnames(n_avg_admix[, 4:ncol(n_avg_admix)]), alpha = 1) +
        scale_fill_manual(values = n_colors)
north
ggsave("lea/distichus/north/K_4_alpha_10_north_subsp.pdf")

#########
## South paleo-island subset plus A. d. dominicensis I and IV: `west` dataset
# south_snmf <- snmf("distichus-south-island-no-miss-50kbLD.geno", K = 1:7, project = "new", repetitions = 50, tolerance = 0.00001, seed = 16, alpha = 10, entropy = TRUE, ploidy = 2)
south_snmf <- load.snmfProject("lea/distichus/south/distichus-south-island-no-miss-50kbLD.snmfProject")
    pdf("lea/distichus/south/snmf_alpha_10_K_1-7.pdf")
        plot(south_snmf, cex = 1.2, col = "lightblue", pch = 19)
    dev.off()
    for (i in 1:7){
        print(paste("For K set to", i, "mean cross entropy is", mean(cross.entropy(south_snmf, K = i))))
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

    s_colors <- c("#003365", "burlywood", "darkorchid2")
    pdf("lea/distichus/south/K3-alpha-10-south-island-subsp-barplot.pdf", width = 8, height = 9)
        barplot(t(s_df[, 6:ncol(s_df)]), col = s_colors, names.arg = s_df$specimen, las = 2, cex.names = 0.3, border = NA)
    dev.off()

    clusters <- grep("V", names(s_df))
    s_avg_admix <- aggregate(s_df[, clusters], by = s_df[, c("locality", "longitude", "latitude")], mean)

basemap + new_scale_fill() +
    geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = s_avg_admix, cols = colnames(s_avg_admix[, 4:ncol(s_avg_admix)]), alpha = 1) +
    scale_fill_manual(values = s_colors)
ggsave("lea/distichus/south/K_3_alpha_10_south_subsp.pdf")

north + new_scale_fill() +
    geom_scatterpie(aes(x = longitude, y = latitude, group = locality), data = s_avg_admix, cols = colnames(s_avg_admix[, 4:ncol(s_avg_admix)]), alpha = 1) +
    scale_fill_manual(values = s_colors)
ggsave("lea/distichus/snmf-pie-plot-hierarchical.pdf")