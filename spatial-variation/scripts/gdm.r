setwd("/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/ibd-ibe/gdm")

# Load R packages
library(tidyverse)
library(adegenet)
library(gdm)

print(sessionInfo())

# Load ddRADseq individual information
columns <- c("specimen", "specimen_pop", "filepath", "subspecies", "population", "complex", "island", "locality", "longitude", "latitude", "note")
RAD_data <- read_table("~/distichus-ddRAD/info/distichus-species-group-Hispaniola-only-popmap-cleaned.tsv", col_names = FALSE)
        colnames(RAD_data) <- columns

# Load linkage-pruned, imputed dataset with ~6,000 SNPs
load("../distichus-filt-imp-genind.RData")
pop(dat_imp) <- RAD_data$locality
# convert to population
pop_imp <- genind2genpop(dat_imp)

## Environmental data
coords <- distinct(RAD_data[, c("longitude", "latitude")])
coords <- as.matrix(coords)

env_df <- read_table("../../selected-environmental-variables.tsv", col_names = TRUE)
env_df <- distinct(env_df)

# Merge environmental data with filtered dataset
# all_data <- data.frame(RAD_data, env_df[, -c(1:3)])

# Choose variables identified by forward selection
# variables <- c("clt_min", "elevation_30sec", "rsds_1981.2010_max_V.2.1", "bio15", "hurs_range", "bio3", "bio9", "bio4", "cmi_range")
variables <- c("clt_min", "elevation_30sec", "bio4", "bio18", "bio9", "cmi_range", "bio3")

# good_env <- all_data[, variables]
good_env <- env_df[, variables]

#################################################
# Run GDM

## Format data for gdm
# sitespecies <- data.frame(all_data$longitude, all_data$latitude, dat_imp@tab)
# sitespecies$site <- rownames(dat_imp@tab)
sitespecies <- data.frame(env_df$longitude, env_df$latitude, pop_imp@tab)
sitespecies$site <- rownames(pop_imp@tab)
pred_env <- good_env
# pred_env$site <- rownames(dat_imp@tab)
pred_env$site <- rownames(pop_imp@tab)
pred_env <- data.frame(pred_env, env_df$longitude, env_df$latitude)
    no_env <- pred_env[, c("site", "env_df.longitude", "env_df.latitude")]

## Create input tables
gdm_tab <- formatsitepair(sitespecies, bioFormat = 1, XColumn = "env_df.longitude", YColumn = "env_df.latitude", predData = pred_env, siteColumn = "site")
gdm_tab_geog <- formatsitepair(sitespecies, bioFormat = 1, XColumn = "env_df.longitude", YColumn = "env_df.latitude", predData = no_env, siteColumn = "site")

## Run GDM

gdm_e <- gdm(gdm_tab, geo = FALSE)
        print(summary(gdm_e))
        save(gdm_e, file = "gdm_env_only.RData")

gdm_eg <- gdm(gdm_tab, geo = TRUE)
        print(summary(gdm_eg))
        save(gdm_eg, file = "gdm_env_geo.RData")
 
gdm_g <- gdm(gdm_tab_geog, geo = TRUE)
        print(summary(gdm_g))
        save(gdm_g, file = "gdm_geo_only.RData")

pdf("gdm_per-loc-full.pdf")
        plot(gdm_eg, plot.layout = c(4, 3))
dev.off()

#gdm_1_imp <- gdm.varImp(gdm_tab, geo = TRUE, nPerm = 999, cores = 12, parallel = T, outFile = "gdm.1.full.imp", fullModelOnly = FALSE) #fullModelOnly argument missing?

# gdm_1_imp_full <- gdm.varImp(gdm_tab, geo = TRUE, nPerm = 999, cores = 12, parallel = T, outFile = "gdm.1.full.imp")
# gdm_1_imp_full_nogeo <- gdm.varImp(gdm_tab, geo = FALSE, nPerm = 999, cores = 12, parallel = TRUE, outFile = "gdm.1.nogeo.imp")
gdm_1_imp_full <- gdm.varImp(gdm_tab, geo = TRUE, nPerm = 999, cores = 12, parallel = T, outFile = "gdm.1.persite.full.imp")
gdm_1_imp_full_nogeo <- gdm.varImp(gdm_tab, geo = FALSE, nPerm = 999, cores = 12, parallel = TRUE, outFile = "gdm.1.persite.nogeo.imp")