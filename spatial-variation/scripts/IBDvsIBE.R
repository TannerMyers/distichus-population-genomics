######################################################################################
# This script includes code for estimating pairwise Fst in the 
# R package `hierfstat`, as well as estimating geographic distance 
# and environmental dissimilarity to quantify the contributions of
# isolation-by-distance and isolation-by-environment to genetic differentiation
# using multivariate statistical approaches.
# Authors: Tanner C. Myers, Pietro L. H. de Mello, Paul M. Hime, and Richard E. Glor
######################################################################################

work_dir <- "/scratch/tcm0036/distichus-ddRAD/analyses/environmental-association/"
setwd(work_dir)

rm(list = ls())

# Load R packages
library(sp)
library(raster)
library(rgeos)
library(rgdal)
library(adegenet)
library(vcfR)
library(hierfstat)
library(adespatial)
library(gdm)
library(vegan)
library(tidyverse)

# ddRADseq individual information
RAD_data <- read_table("/home/tcm0036/distichus-ddRAD/info/distichus-popmap-cleaned-master.tsv", col_names = TRUE)
snmfK9 <- read_table("/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/lea/sNMF_K9_ancestry_coefficients_cluster-specimen-cleaned.tsv")
    names(snmfK9)[names(snmfK9) == "V1"] <- "Sample_ID_pop"
RAD_data2 <- RAD_data %>% filter(Sample != "4481" & Sample != "10085" & Sample != "10086") %>%
    filter(!(Taxon %in% c("websteri", "marron", "caudalis", "altavelensis"))) %>%
    filter(Island %in% "Hispaniola")  %>%
    filter(!grepl("_rep", Sample)) # Temporary, eventually I need to deal with the technical replicates in the dataset
RAD_data2 <- inner_join(snmfK9, RAD_data2)
RAD_data2 <- RAD_data2[RAD_data2$cluster != 6 & RAD_data2$cluster != 8, ]

# Estimate distance matrices
## Geographic distance
coords <- distinct(RAD_data2[, c("Longitude", "Latitude")])
coords <- as.matrix(coords)

geog_dist <- sp::spDists(x = coords, y = coords, longlat = TRUE)

## Environmental distance
env_df <- read_csv("mmrr/raw-environmental-variables.csv", col_names = TRUE)

env_df <- inner_join(as_tibble(coords), env_df)

# Now, do the same thing, but for raw environmental variables
env_dist <- as.matrix(dist(env_df[, 3:ncol(env_df)]))

all_data <- inner_join(RAD_data2, env_df)

latlong <- SpatialPoints(RAD_data2[, c("Longitude", "Latitude")], proj4string = CRS("+proj=longlat"))

distichus_xy <- as.data.frame(spTransform(latlong, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")))

dbmem_all <- dbmem(distichus_xy)
    load("rda/dbMEM_distichus-subspecies.RData")
dbmem_keep <- dbmem_all[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 20, 21, 24)]

vcf_dir <- "/scratch/tcm0036/distichus-ddRAD/analyses/environmental-association/vcfs-by-chrom/"

# Loop over individual chromosome vcf files to
## (1) Perform redundancy analysis
## (2) Perform generalized dissimilarity modeling analysis
for (i in list.files(vcf_dir, pattern = ".vcf")){

setwd(work_dir)

# Load relevant vcf
vcfR <- read.vcfR(paste0(vcf_dir, i))
# convert vcfR object to genind object
dat <- vcfR2genind(vcfR)
dat2 <- dat[indNames(dat) %in% RAD_data2$Sample_ID_pop]
pop(dat2) <- RAD_data2$Locality
    assign(paste0(str_replace(i, ".vcf", ""), "_dat"), dat2)

gen <- vcfR2genlight(vcfR)
gen2 <- gen[gen@ind.names %in% RAD_data2$Sample_ID_pop]
pop(gen2) <- RAD_data2$Locality
    assign(paste0(str_replace(i, ".vcf", ""), "_gen"), gen2)

    # data <- genind2hierfstat(dat = dat2, pop = dat2@pop)
    # data <- as.data.frame(cbind(dat2@pop, data[, c(2:ncol(data))]))

    # gen_diff <- pairwise.WCfst(dat = data, diploid = TRUE)
    # #assign(gen_diff, str_replace(i, ".vcf", "_fst"))
    #     save(gen_diff, file = paste0("vcfs-by-chrom/", str_replace(i, ".vcf", ""), "-Fst.RData"))
    
    # gen_dist <- gen_diff/(1 - gen_diff)
    # class(gen_dist)
    # isSymmetric(gen_dist) # Needs to be TRUE

# rda_distichus_1 <- rda(dat2@tab ~., data = all_data[, 22:ncol(all_data)], scale = TRUE)
#     save(rda_distichus_1, file = paste0("rda/", str_replace(i, ".vcf", ""), "_rda_1.RData"))
# rda_distichus_1r2 <- RsquareAdj(rda_distichus_1)

# levels(all_data$cluster) <- c("ignigularis", "east dominicensis", "South paleo-island", "properus", "west dominicensis", "ravitergum")
# clust <- all_data$cluster
# bg <- c("#cf5d34", "#634102", "#993cc8", "gray80", "#c6b89c", "#fddc1f")

# pdf(paste0("rda/", str_replace(i, ".vcf", ""), "-full-plot-simple.pdf"))
#     plot(rda_distichus_1, type = "n", scaling = 3)
#     #points(rda_distichus_1, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3, choices = c(1, 2)) # SNPs in this case
#     points(rda_distichus_1, display = "sites", pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = bg[as.factor(clust)], choices = c(1, 2)) # lizards
#     #text(rda_distichus_1, scaling = 3, display = "bp", col = "#0868ac", cex = 1) # predictor variables
#     legend("topright", legend = levels(clust), bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bg)
# dev.off()

# # Forward selection of environmental variables
# rda_distichus_fwd <- forward.sel(dat2@tab, all_data[, 22:ncol(all_data)], adjR2thresh = rda_distichus_1r2, nperm = 9999, alpha = 0.01)
#     save(rda_distichus_fwd, file = paste0("rda/", str_replace(i, ".vcf", ""), "_fwdsel.RData"))
    load(paste0("rda/", str_replace(i, ".vcf", ""), "_fwdsel.RData"))

variables <- rda_distichus_fwd$variables
good_env <- all_data[, variables]

# rda_distichus_2 <- rda(dat2@tab ~., data = good_env, scale = TRUE)
#     assign(paste0(str_replace(i, ".vcf", ""), "_rda_2"), rda_distichus_2)

# pdf(paste0("rda/", str_replace(i, ".vcf", ""), "fwd-sel-plot.pdf"))
#     plot(rda_distichus_2, type = "n", scaling = 3)
#     points(rda_distichus_2, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3, choices = c(1, 2)) # SNPs in this case
#     points(rda_distichus_2, display = "sites", pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = bg[as.factor(clust)], choices = c(1, 2)) # lizards
#     text(rda_distichus_2, scaling = 3, display = "bp", col = "#0868ac", cex = 1) # predictor variables
#     legend("bottomright", legend = levels(clust), bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bg)
# dev.off()

# ######## Variance partitioning

# vp <- varpart(dat2@tab, good_env, dbmem.1$dbMEM_red_model)
# png(paste0("rda/", str_replace(i, ".vcf", ""), "venn_partitions.png"), width = 6, height = 6, units = "in", res = 600)
#     plot(vp, bg = c("red", "blue"), digits = 1, Xnames = c('Env', 'Space'))
# dev.off()

# ######## Constrained on Space
# rda_distichus_cond <- rda(dat2@tab, good_env, dbmem.1$dbMEM_red_model)
#    save(rda_distichus_cond, file = paste0("rda/", str_replace(i, ".vcf", ""), "_rda_cond.RData"))

# null_model <- rda(dat2@tab ~ Condition(as.matrix(dbmem.1$dbMEM_red_model)), data = all_data[, 22:ncol(all_data)])

# full_model <- rda(dat2@tab ~ . + Condition(as.matrix(dbmem.1$dbMEM_red_model)), data = all_data[, 22:ncol(all_data)])

# ordi <- ordistep(null_model, scope = formula(full_model), direction = "forward")
#     save(ordi, file = paste0("rda/", str_replace(i, ".vcf", ""), "constrained-variable-selection.RData"))

# Run GDM


setwd("gdm")
sitespecies <- data.frame(all_data$Longitude, all_data$Latitude, dat2@tab)
sitespecies$site <- rownames(dat2@tab)
pred_env <- good_env
pred_env$site <- rownames(dat2@tab)
pred_env <- data.frame(pred_env, all_data$Longitude, all_data$Latitude)
    no_env <- pred_env[, c(7, 8, 9)]

# Format data for GDM
gdm_tab <- formatsitepair(sitespecies, bioFormat = 1, XColumn = "all_data.Longitude", YColumn = "all_data.Latitude", predData = pred_env, siteColumn = "site")
gdm_tab_geog <- formatsitepair(sitespecies, bioFormat = 1, XColumn = "all_data.Longitude", YColumn = "all_data.Latitude", predData = no_env, siteColumn = "site")

# Run GDM
gdm_e <- gdm(gdm_tab, geo = FALSE)
gdm_eg <- gdm(gdm_tab, geo = TRUE)
gdm_g <- gdm(gdm_tab_geog, geo = TRUE)

# Run GDM with variable selection
gdm_1_imp_full <- gdm.varImp(gdm_tab, geo = TRUE, nPerm = 999, cores = 12, parallel = TRUE, outFile = paste0(str_replace(i, ".vcf", ""), "gdm.1.full.imp")) # fullModelOnly argument not available?
gdm_1_imp_full_nogeo <- gdm.varImp(gdm_tab, geo = FALSE, nPerm = 999, cores = 12, parallel = TRUE, outFile = paste0(str_replace(i, ".vcf", ""), "gdm.1.nogeo.imp"))

pred_env2 <- subset(pred_env, select = -c(CHELSA_bio10_02, CHELSA_bio10_03, CHELSA_bio10_04, CHELSA_bio10_07))
gdm_tab2 <- formatsitepair(sitespecies, bioFormat=1, XColumn = "all_data.Longitude", YColumn = "all_data.Latitude", predData = pred_env2, siteColumn = "site")

gdm(gdm_tab2, geo = TRUE)
}