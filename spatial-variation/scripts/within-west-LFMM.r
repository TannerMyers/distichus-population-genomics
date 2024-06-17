setwd("/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/gea/lfmm/LD-pruned/west")

library(tidyverse)
library(LEA)
library(tess3r)
library(terra)
library(qvalue)
library(lfmm)

print(sessionInfo())

## Load specimen information
west_sub <- read_table("~/distichus-ddRAD/info/south-island-dataset-master.tsv",
	col_names = c('specimen', 'specimen_pop', 'filepath', 'subspecies', 'population', 'complex', 'island', 'locality', 'longitude', 'latitude', 'note'))

vcf <- "west_subset_imp_LD.recode.vcf"

## Prepare data for LFMM: first, SNPs
# vcf2lfmm(vcf, output.file = "west_subset_imp_LD.recode.lfmm") # also makes .geno file
lfmm <- read.lfmm("west_subset_imp_LD.recode.lfmm")
	print(dim(lfmm))

## Prepare data for LFMM: now environmental variables
env <- read_table("../../../fwd-select-env-pc1-3.tsv")
west_env <- env[env$specimen_pop %in% west_sub$specimen_pop,]
        dim(west_env)
        sum(is.na(west_env))
west_env <- west_env[, c("Axis1", "Axis2", "Axis3")]
        
write.env(west_env, "west-lfmm_pca_env.env")
env <- read.env(input.file = "west-lfmm_pca_env.env")

lfmm_r <- lfmm_ridge(Y = lfmm, X = env, K = 4)
	save(lfmm_r, file = "west-lfmm-ridge.RData")

K3_pv <- lfmm_test(Y = lfmm, X = env, lfmm = lfmm_r, calibrate = "gif")
	save(K3_pv, file = "west-lfmm-test.RData")
	print(K3_pv)
	print(names(K3_pv))
	print(K3_pv$gif)
