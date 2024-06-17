
setwd("/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/gea/lfmm/LD-pruned")

library(tidyverse)
library(LEA)
library(tess3r)
library(terra)
library(qvalue)
library(lfmm)

print(sessionInfo())

#########################

## imputed VCF (139,964 SNPs for 191 specimens)
# vcf <- "/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/vcf/distichus-subsp-only/imputation/distichus_geno0.8_ind0.8_maf0.05_dp4_imputed.vcf"
#vcf <- "/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/vcf/distichus-subsp-only/imputation/distichus_imp_LD.vcf"

## Load specimen information
columns <- c("specimen", "specimen_pop", "filepath", "subspecies", "population", "complex", "island", "locality", "longitude", "latitude", "note")
RAD_data <- read_table("~/distichus-ddRAD/info/distichus-species-group-Hispaniola-only-popmap-cleaned.tsv", col_names = FALSE)
         colnames(RAD_data) <- columns

## Prepare data for LFMM: First, SNPs
#vcf2lfmm(vcf, output.file = "distichus_imp_LD.lfmm") # also makes .geno file
lfmm <- read.lfmm("distichus_imp_LD.lfmm")
write.env(env, "lfmm_pca_env.env")
env <- read.env(input.file = "lfmm_pca_env.env")

######## Running lfmm 2.0
## Run lfmm using least squares estimation with a ridge penalty
lfmm_pca_r <- lfmm_ridge(Y = lfmm, X = env, K = 5)
    save(lfmm_pca_r, file = "LFMM2/lfmm-ridge-3pcs.RData")
    load("lfmm/LD-pruned/lfmm-ridge-3pcs.RData") # local

# Perform association test between each column of response matrix Y
# and explanatory matrix X, while correcting for unobserved confounders 
# (hence the latent)
K5_pca.pv <- lfmm_test(Y = lfmm, X = env, lfmm = lfmm_pca_r, calibrate = "gif")
	save(K5_pca.pv, file = "LFMM2/lfmm-pca-test.RData")
    load("lfmm/LD-pruned/lfmm-pca-test.RData")
print(names(K5_pca.pv))
# [1] "B"                 "epsilon.sigma2"    "B.sigma2"         
# [4] "score"             "pvalue"            "gif"              
# [7] "calibrated.score2" "calibrated.pvalue"

## Look at genomic inflation factors; should be ~1 if appropriately calibrated
print(K5_pca.pv$gif)
#       V1       V2       V3 
# 2.485622 2.002009 1.795677 

for (pc in 1:3){
    png(paste0("P-values-PC", pc, ".png"))
        par(mfrow = c(2, 1))
        hist(K5_pca.pv$pvalue[, pc], main = "Unadjusted p-values")
        hist(K5_pca.pv$calibrated.pvalue[, pc], main = "GIF-adjusted p-values")
    dev.off()
}

## Manually adjust GIF to achieve appropriate p-value distribution
## (i.e., flat with peak around 0)
## PC1
zscore_1 <- K5_pca.pv$score[, 1]
gif_1 <- K5_pca.pv$gif[1] 

new_gif1 <- 2.6

adj_pv1 <- pchisq(zscore_1^2/new_gif1, df=1, lower = FALSE)

pdf("p-values-PC1-distributions.pdf", width = 12, height = 8)
par(mfrow=c(1, 3))
hist(K5_pca.pv$pvalue[, 1], main="Unadjusted p-values")        
hist(K5_pca.pv$calibrated.pvalue[, 1], main=paste("GIF-adjusted p-values (GIF =", gif_1, ")"))
hist(adj_pv1, main=paste("REadjusted p-values (GIF=", new_gif1, ")"))
dev.off()

## Manually adjust GIF to achieve appropriate p-value distribution
## (i.e., flat with peak around 0)
## PC2
zscore_2 <- K5_pca.pv$score[, 2]
gif_2 <- K5_pca.pv$gif[2] 

new_gif2 <- 2.05

adj_pv2 <- pchisq(zscore_2^2/new_gif2, df=1, lower = FALSE)

pdf("p-values-PC2-distributions.pdf", width = 12, height = 8)
par(mfrow=c(1, 3))
hist(K5_pca.pv$pvalue[, 2], main="Unadjusted p-values")        
hist(K5_pca.pv$calibrated.pvalue[, 2], main=paste("GIF-adjusted p-values (GIF =", gif_2, ")"))
hist(adj_pv2, main=paste("REadjusted p-values (GIF=", new_gif2, ")"))
dev.off()

## Manually adjust GIF to achieve appropriate p-value distribution
## (i.e., flat with peak around 0)
## PC3
zscore_3 <- K5_pca.pv$score[, 3]
gif_3 <- K5_pca.pv$gif[3] 

new_gif3 <- 1.85

adj_pv3 <- pchisq(zscore_3^2/new_gif3, df=1, lower = FALSE)

pdf("p-values-PC3-distributions.pdf", width = 12, height = 8)
par(mfrow=c(1, 3))
hist(K5_pca.pv$pvalue[, 3], main="Unadjusted p-values")        
hist(K5_pca.pv$calibrated.pvalue[, 3], main=paste("GIF-adjusted p-values (GIF =", gif_3, ")"))
hist(adj_pv3, main=paste("REadjusted p-values (GIF=", new_gif3, ")"))
dev.off()

## Now convert p to q values and implement FDR correction

## First, for the calibrated p-values
K5_qv <- qvalue(K5_pca.pv$calibrated.pvalue)$qvalues
    length(which(K5_qv < 0.05)) # how many SNPs have a false discovery rate less than 5%? 

## Now, for my adjusted values
qv_1 <- qvalue(adj_pv1)$qvalues
    length(which(qv_1 < 0.05))
qv_2 <- qvalue(adj_pv2)$qvalues
    length(which(qv_2 < 0.05))
qv_3 <- qvalue(adj_pv3)$qvalues
    length(which(qv_3 < 0.05))

## Both of these are within a few dozen of each other, between 520 and 550, as a result, I will go with the non-adjusted calibrated p-values 

vcfsnp <- read_table(col_names = FALSE, "distichus_imp_LD.vcfsnp")
snp_names <- paste0(vcfsnp$X1, "_", vcfsnp$X2)
marker_pos <- paste0(snp_names, "_", 1:39427)

#all <- cbind.data.frame(snp_names, qv_1, qv_2, qv_3, 1:39427)
all <- cbind.data.frame(snp_names, marker_pos, K5_qv, 1:39427)
colnames(all) <- c("snp", "marker_pos", "pc1", "pc2", "pc3", "integer")

## subset to just candidates/outliers
all_cand <- all[all$pc1 < 0.05 |
                    all$pc2 < 0.05 |
                    all$pc3 < 0.05,]
dim(all_cand)
all_cand <- na.omit(all_cand)
write.csv(all_cand, "lfmm-candidate-SNPs.csv", row.names = FALSE, quote = FALSE)

###################################################################################################
## Now, run LFMM with original MCMC recipe
lfmm <- lfmm("distichus_imp_LD.lfmm", "lfmm_pca_env.env", K = 5, project = "new", repetitions = 5)

#dist_lfmm <- load.lfmmProject("distichus_imp_LD_lfmm_env.lfmmProject")
raw_lfmm <- load.lfmmProject("distichus_imp_LD_lfmm_env.lfmmProject")
pca_lfmm <- load.lfmmProject("distichus_imp_LD_lfmm_pca_env.lfmmProject")

## Calculate adjusted p-values and determine outliers from forward-selected environmental variables
for (p in 1:7){
    pvals <- lfmm.pvalues(raw_lfmm, K = 5, d = p)
    zscores <- z.scores(raw_lfmm, K = 5, d = p)
    zscores_med <- apply(zscores, 1, median, na.rm = TRUE)

    new_gif <- 1
    adj_pval <- pchisq(zscores_med^2/new_gif, df = 1, lower = FALSE)

    ## Select outliers based on q-values
    p_adjusted <- p.adjust(adj_pval, method = "fdr", n = 39427)
    alpha = 0.05        ## false discovery rate
    outliers <- which(p_adjusted < alpha)
    length(outliers)

    ggplot() + geom_histogram(aes(adj_pval)) + labs(x = "p-values", y = "Count") +
        ggtitle(paste0("P-values -- ", "Variable-", p)) + theme_minimal()
    ggsave(paste0("Variable-", p, " p-values.png"))
    ggplot() + geom_histogram(aes(p_adjusted)) + labs(x = "Adj. p-values", y = "Count") +
        ggtitle(paste0("Adjusted p-values -- ", "Variable-", p)) + theme_minimal()
    ggsave(paste0("Variable-", p, " p-values.png"))

    ## Keep objects in environment
    assign(paste0("pvals_var", p), pvals)
    assign(paste0("adjpval_var", p), p_adjusted)
    assign(paste0("outliers_var", p), outliers)
}

# PCA p-values
for (p in 1:3){
    pvals <- lfmm.pvalues(pc_lfmm, K = 5, d = p)
    zscores <- z.scores(pc_lfmm, K = 5, d = p)
    zscores_med <- apply(zscores, 1, median, na.rm = TRUE)

    new_gif <- 1
    adj_pval <- pchisq(zscores_med^2/new_gif, df = 1, lower = FALSE)

    ## Select outliers based on q-values
    p_adjusted <- p.adjust(adj_pval, method = "fdr", n = 39427)
    alpha <- 0.05        ## false discovery rate
    outliers <- which(p_adjusted < alpha)
    length(outliers)

    ggplot() + geom_histogram(aes(adj_pval)) + labs(x = "p-values", y = "Count") +
        ggtitle(paste0("P-values -- ", "Principal Component -", p)) + theme_minimal()
    ggsave(paste0("PC-", p, "-p-values.png"))
    ggplot() + geom_histogram(aes(p_adjusted)) + labs(x = "Adj. p-values", y = "Count") +
        ggtitle(paste0("Adjusted p-values -- ", "Principal Component -", p)) + theme_minimal()
    ggsave(paste0("PC-", p, "-adj-p-values.png"))

    ## Keep objects in environment
    assign(paste0("pvals_var", p), pvals)
    assign(paste0("adjpval_var", p), p_adjusted)
    assign(paste0("outliers_var", p), outliers)
}

vcfsnp <- read_table(col_names = FALSE, "distichus_imp_LD.vcfsnp")
snp_names <- paste0(vcfsnp$X1, "_", vcfsnp$X2)

# pc1 <- snp_names[outliers_var1]
# pc2 <- snp_names[outliers_var2]
# pc3 <- snp_names[outliers_var3]

all_pvals <- cbind.data.frame(snp_names, adjpval_var1, adjpval_var2, adjpval_var3, 1:39427)
colnames(all_pvals) <- c("snp", "pc1", "pc2", "pc3", "integer")

## subset to just candidates/outliers
all_cand <- all_pvals[all_pvals$pc1 < 0.05 |
                      all_pvals$pc2 < 0.05 |
                      all_pvals$pc2 < 0.05,]
write.csv(all_cand, "lfmm-candidate-SNPs.csv")