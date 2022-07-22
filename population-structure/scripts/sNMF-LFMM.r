####################################################################################
# This script includes code for estimating population structure with
# sparse non-negative matrix factorization (sNMF) in the R package `LEA`
# and then using the number of clusters to inform genotype-environment
# association analysis with the univariate approach latent factor mixture
# models implemented in the R package `LFMM`.
# Authors: Tanner C. Myers, Pietro L. H. de Mello, Paul M. Hime, and Richard E. Glor
####################################################################################

rm(list = ls())

library(tidyverse)
library(LEA)
library(tess3r)
library(terra)
library(qvalue)
library(lfmm)

#########################

setwd("/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/lea")

# Load specimen information
RAD_data <- read_table("/home/tcm0036/distichus-ddRAD/info/distichus-popmap-cleaned-master.tsv",
                     col_names = TRUE)

RAD_data2 <- RAD_data %>% filter(Sample != "4481" & Sample != "10085" & Sample != "10086") %>% 
    filter(!(Taxon %in% c("websteri", "marron", "caudalis", "altavelensis"))) %>% 
    filter(Island %in% "Hispaniola") %>%
    filter(!grepl("_rep", Sample)) # Temporary, eventually I need to deal with the technical replicates in the dataset

################################################################
## Infer ancestry proportions with sNMF

vcf <- "../distichus-spgrp-no-missing-LD-pruned-informative-names.vcf"
vcf2geno(vcf, output.file = "distichus-spgrp-no-missing-LD-pruned.geno")
geno <- read.geno("distichus-spgrp-no-missing-LD-pruned.geno")
    dim(geno)

obj.snmf <- snmf("distichus-spgrp-no-missing-LD-pruned.geno", K = 1:10, project = "new",
     repetitions = 50, tolerance = 0.00001, seed = 16,
     alpha = 100, entropy = TRUE, ploidy = 2)
     obj.snmf <- load.snmfProject("distichus-spgrp-no-missing-LD-pruned.snmfProject")

pdf("snmf_K1-10_50-replicates.pdf")
    plot(obj.snmf, cex = 1.2, col = "lightblue", pch = 19)
dev.off()

for (i in 1:10){
     print(mean(cross.entropy(obj.snmf, K = i)))
}

# For optimal value of K, plot ancestry coefficients
for (i in 9:10){
    k = i
    ce <- cross.entropy(obj.snmf, K = i) # K = 10 determined optimal
    best_run <- which.min(ce)

    qmatrix <- LEA::Q(obj.snmf, run = best_run, K = k)
    Qmatrix  <- tess3r::as.qmatrix(qmatrix)

    # get cluster assignments per individual
    cluster <- apply(qmatrix, 1, which.max)

    sNMF_ancestry_coefficients_cluster <- as.data.frame(cbind(RAD_data$Sample_ID_pop, cluster, Qmatrix))
        write_delim(x = sNMF_ancestry_coefficients_cluster,
            file = paste0("sNMF_K", k, "_ancestry_coefficients_cluster-specimen-cleaned.tsv"),
            col_names = TRUE, delim = "\t")

    snmf <- sNMF_ancestry_coefficients_cluster[order(as.numeric(sNMF_ancestry_coefficients_cluster$cluster)), ]

    if (k == 10){
        pdf(paste0("sNMF-K", k, "-barplot-specimen-cleaned.pdf"))
            barplot(t(snmf[,3:ncol(snmf)]),
                col = c("black", "#c6b89c", "#1a1ac2", "gray80", "#634102",
                        "#993cc8", "black", "#04e779", "#cf5d34", "#fddc1f"),
                names.arg = snmf$V1, las = 2, cex.names = 0.3, border = NA)
        dev.off()
    } else {
        pdf(paste0("sNMF-K", k, "-barplot-specimen-cleaned.pdf"))
            barplot(t(snmf[,3:ncol(snmf)]),
                col = c("gray80", "#634102", "#c6b89c", "#04e779",
                        "#fddc1f", "#993cc8", "#cf5d34", "black", "#1a1ac2"),
                names.arg = snmf$V1, las = 2, cex.names = 0.3, border = NA)
        dev.off()

        #pdf(paste0("sNMF-K", k, "-spatial.pdf"))
        #    plot(x = , RAD_data[,c('Longitude', 'Latitude')], 
        #        col.palette = c("gray80", "#634102", "#c6b89c", "#04e779",
        #                     "#fddc1f", "#993cc8", "#cf5d34", "black", "#1a1ac2"),
        #     interpolation.model = FieldsKrigModel(10))
        # dev.off()
    }
}

################################################################
# Now to look at environment-genotype associations with latent
# factor mixture models (LFMM)

# First, need to generate a file in the lfmm format
# geno2lfmm(input.file = "distichus-spgrp-no-missing-LD-pruned.geno", output.file = "distichus-spgrp-no-missing-LD-pruned.lfmm")
# lfmm <- read.lfmm("distichus-spgrp-no-missing-LD-pruned.lfmm")
#     dim(lfmm)

# The object `lfmm` contains more samples that I am interested in performing
# environmental analyses on at the moment. To rectify this, I will try to
# index the matrix to include only rows corresponding to specimens
# in the RAD_data2 dataframe
# matched <- RAD_data$Sample %in% RAD_data2$Sample
#     lfmm2 <- lfmm[matched, ]
#     write.lfmm(R = lfmm2, output.file = "distichus-spgrp-Hispaniola-only.lfmm")
lfmm <- read.lfmm("distichus-spgrp-Hispaniola-only.lfmm")

# Now, we need an `env` object that contains values for environmental variables for each locality
# clim <- raster::stack(list.files(path = "~/distichus-spatial-variation/data/chelsa_new/", 
#                                 pattern = ".asc", full.names = TRUE))

# fill.na <- function(x) {
#   center = 0.5 + (width*width/2)
#   if(is.na(x)[center]) {
#     return(mean(x, na.rm = TRUE))
#   } else {
#     return(x[center])
#   }
# }

# width <- 5
# # for (i in c(1:15)){
# #     var <- terra::focal(clim[[i]], w = matrix(1, width, width), fun = fill.na, na.rm = FALSE, pad = TRUE)
# #     assign(paste0("var", i), var)
# # }

# for (i in c(1:15)){
#     var <- terra::focal(clim[[i]], w = matrix(1, width, width), 
#     fun = fill.na, na.rm = FALSE, pad = TRUE)

#     ## Because bio8 and bio9 have been removed,
#     ## fix names when assigning variables
#     if (i >= 8){
#         j <- i + 2
#         assign(paste0("bio", j), var)
#     } else {
#         assign(paste0("bio", i), var)
#     }
# }

# clim2 <- stack(bio1, bio2, bio3, bio4, bio5,
#                 bio6, bio7, bio10, bio11, bio12,
#                 bio13, bio14, bio15, bio16, bio17)
#     names(clim2) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6",
#                 "bio7", "bio10", "bio11", "bio12", "bio13", "bio14",
#                 "bio15", "bio16", "bio17")

# clim_vals <- terra::extract(clim2, RAD_data2[, c("Longitude", "Latitude")])
#     write.env(R = clim_vals, output.file = "CHELSA-Hispaniola-distichus-only.env")
env <- read.env("CHELSA-Hispaniola-distichus-only.env")

# Run lfmm
## Instead of using the MCMC approach of the initial LFMM introduced by
## Frichot et al. (2013), I am using the least squares regression approach
## with a ridge penalty introduced by Caye et al. (2019).
mod.K9.lfmm <- lfmm_ridge(Y = lfmm, X = env, K = 9)

# Perform association test between each column of response matrix Y
# and explanatory matrix X, while correcting for unobserved confounders 
# (hence the latent)
K9.pv <- lfmm_test(Y = lfmm, X = env, lfmm = mod.K9.lfmm, calibrate = "gif")
    names(K9.pv)

## Look at genomic inflation factors; should be ~1 if appropriately calibrated
K9.pv$gif

## Look at p-value distributions
pdf("LFMM-K9-pvalue-distribution.pdf")
    par(mfrow=c(2,1))
    hist(K9.pv$pvalue[, 1], main = "Unadjusted p-values")
    hist(K9.pv$calibrated.pvalue[, 1], main = "GIF-adjusted p-values")
dev.off()

## Convert p-values to q-values that provide measure of each SNP's significance
K9.qv <- qvalue(K9.pv$calibrated.pvalue)$qvalues
    ## How many SNPs have an FDR < 10%?
    length(which(K9.qv < 0.1))

# Find number of correlated SNPs for each bioclimatic variable

theta <- 4.5
pdf("LFMM-K9-number-of-SNPs-per-variable.pdf")
    barplot(apply(K9.pv$score, 2, FUN = function(x) sum(abs(na.omit(x)) > theta)),
            ylab = "Number of hits", las = 3,
            names.arg = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7",
                        "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17")
            )
dev.off()

## Importance of each bioclimatic variable
importance <- apply((abs(na.omit(K9.pv$score)) > theta) * na.omit(K9.pv$score), 2, FUN = function(x) mean(x^2))
pdf("LFMM-K9-variable-importance.pdf")
    barplot(importance, col = "blue4", xlab = "Bioclimatic variables", ylab = "Importance", las = 3,
            names.arg = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7",
                        "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17")
            )
dev.off()

hit <- (abs(na.omit(K9.pv$score)) > theta)
result <- NULL
for (i in 1:15){
    if (sum(hit[, i]) > 0){
        y <- env[, i]
        data <- as.data.frame(lfmm[, hit[, i]])
        mod <- lm(y ~ ., data = data)
        smr <- summary(mod)
        result <- rbind(result, c(i, smr$r.squared, smr$fstatistic, pf(smr$fstatistic[1], smr$fstatistic[2], smr$fstatistic[3], low = FALSE)))
    } else {
        result <- rbind(result, c(i, 0, 0, 0, 0, 1))
    }
}
rownames(result) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17")

## Multiple R squared statistic for each bioclimatic variable
pdf("LFMM-K9-multiple-Rsquared-per-variable.pdf")
    barplot(result[,2], ylab = "Multiple R-squared", las = 3)
dev.off()

result <- cbind(result[,-1], apply(hit, 2, sum))

colnames(result) <- c("MultipleRsquared", "F", "df1", "df2", "pvalue", "hits")
result.df <- data.frame(round(result, digit = 4))