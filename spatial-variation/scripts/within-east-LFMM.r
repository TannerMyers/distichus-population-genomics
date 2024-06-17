setwd("/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/gea/lfmm/LD-pruned/east")

library(tidyverse)
library(LEA)
library(tess3r)
library(terra)
library(qvalue)
library(lfmm)

print(sessionInfo())

## Load specimen information
east_sub <- read_table("~/distichus-ddRAD/info/north-island-dataset-master.tsv",
	col_names = c('specimen', 'specimen_pop', 'filepath', 'subspecies', 'population', 'complex', 'island', 'locality', 'longitude', 'latitude', 'note'))

vcf <- "east_subset_imp_LD.recode.vcf"

## Prepare data for LFMM: first, SNPs
# vcf2lfmm(vcf, output.file = "east_subset_imp_LD.recode.lfmm") # also makes .geno file
lfmm <- read.lfmm("east_subset_imp_LD.recode.lfmm")
	print(dim(lfmm))

## Prepare data for LFMM: now environmental variables
env <- read_table("../../../fwd-select-env-pc1-3.tsv")
east_env <- env[env$specimen_pop %in% east_sub$specimen_pop,]
        dim(east_env)
        sum(is.na(east_env))
east_env <- east_env[, c("Axis1", "Axis2", "Axis3")]
        
write.env(east_env, "east-lfmm_pca_env.env")
env <- read.env(input.file = "east-lfmm_pca_env.env")
        
write.env(env, "east-lfmm_pca_env.env")
env <- read.env(input.file = "east-lfmm_pca_env.env")

lfmm_r <- lfmm_ridge(Y = lfmm, X = env, K = 4)
	save(lfmm_r, file = "east-lfmm-ridge.RData")

K4_pv <- lfmm_test(Y = lfmm, X = env, lfmm = lfmm_r, calibrate = "gif")
	save(K4_pv, file = "east-lfmm-test.RData")
	print(K4_pv)
	print(names(K4_pv))
	print(K4_pv$gif)
