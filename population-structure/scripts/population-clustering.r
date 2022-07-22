setwd("/mmfs1/scratch/tcm0036/distichus-ddRAD/analyses/population-structure")

rm(list = ls())

library(tidyverse)
library(RColorBrewer)
library(adegenet)
library(mclust)
library(MASS)
library(ggforce)
library(cluster)
library(factoextra)
library(gridExtra)
library(randomForest)
library(PCDimension)
library(tsne)
library(vcfR)
library(LEA)
library(tess3r)

# Load specimen data, VCF, and do sub-setting/format conversion
RAD_data <- read_table("/home/tcm0036/distichus-ddRAD/info/distichus-popmap-cleaned-master.tsv", col_names = TRUE)
RAD_data2 <- RAD_data %>% filter(Sample != "4481" & Sample != "10085" & Sample != "10086") %>% 
    filter(!grepl("_rep", Sample))
    ## Take specimens that are removed after sub-setting
    test <- RAD_data$Sample_ID_pop %in% RAD_data2$Sample_ID_pop
        test <- c(TRUE, test)

## Load VCF into R
vcfR <- read.vcfR("distichus-spgrp-no-missing-LD-pruned-informative-names.vcf")
    ## Subset VCF based on above sub-setting
    vcfR <- vcfR[, test]

    ## Convert VCF to genind object
    data <- vcfR2genind(vcfR)
        pop(data) <- RAD_data2$PopID
    ## Convert VCF to genlight object
    gen <- vcfR2genlight(vcfR)
        pop(gen) <- RAD_data2$PopID

    ## Scale genind object
    data_scaled <- scaleGen(data, center = FALSE, scale = FALSE, NA.method = c("mean"), nf)
    data_scaled <- scaleGen(data, center = FALSE, scale = FALSE)

################################################################
# Perform Principal Component Analysis

## This performs PCA on the genlight object `gen` and requires manually choosing the number of principal components to retain
gen_pca <- glPca(gen) # I chose 6
gen_pca_scores <- as_tibble(as.data.frame(gen_pca$scores))
gen_pca_scores$pop <- gen@pop
gen_pca_scores$ind <- gen@ind.names

    ## Set color palette
    colorCount <- length(unique(gen@pop))
    myColors <- colorRampPalette(brewer.pal(9, "Set1"))(colorCount)
    
    p <- ggplot(gen_pca_scores, aes(x = PC1, y = PC2, color = pop))
        p <- p + geom_point(size = 2) + scale_colour_manual(values = myColors)
        p <- p + geom_hline(yintercept = 0)
        p <- p + geom_vline(xintercept = 0)
        p <- p + theme_bw()
    p

    # This performs PCA on the genind object `data`

# Perform Discriminant Analysis of Principal Components

for (i in 1:50){
    set.seed(i)
    grp <- find.clusters(gen, max.n.clust = 20, n.pca = 300, choose.n.clust = FALSE, criterion = "goodfit", #"diffNgroup", 
                        stat = "BIC")
    print(grp$size)
}
## The most often recovered number of clusters is 7 with 27 out of 50 seeds yielding 7 clusters. 
## 6 was the second most common with 13, 8 had 4, 9 had 3, and, oddly, 16 had 2
#  42 34 19 47 24 21 50 30
#  21 42 24 30 47 20 83
#  50 33 47 21 24 42 50
#  83 20 21 27 24 50 42
#  47 20 42 30 21 24 83
#  42 50 47 24 50 33 21
#  19 21 50 21 42 47 12 31 24
#  42 24 83 50 47 21
#  24 50 50 21 47 42 33
#  33 21 24 47 42 50 50
#  47 24 21 42 50 83
#  50 30 47 34 24 19 42 21
#  47 30 24 34 50 19 42 21
#  20 21 42 24 30 30 19 47 34
#  50 33 42 68 50 24
#   20 15  4 12 19  5  7 30 20 34  3 10 30 21 10 27
#    8 10 22 11 12  9 10 11 11 19 15 47 31 10 20 21
#  83 24 50 21 47 42
#  83 47 21 50 42 24
#  33 50 47 50 24 42 21
#  33 42 50 21 74 47
#  21 24 42 34 50 30 47 19
#  50 33 42 47 24 21 50
#  21 47 33 24 42 50 50
#  33 24 21 47 42 50 50
#  47 21 83 10 50 24 32
#  21 50 24 42 33 50 47
#  42 47 83 24 21 50
#  42 50 21 47 24 83
#  21 42 50 33 50 47 24
#  24 50 42 97 21 33
#  21 39 24 83 42 47 11
#  21 24 42 50 50 47 33
#  83 21 47 24 50 42
#  47 24 50 42 33 50 21
#  47 42 83 50 24 21
#  24 68 50 50 33 42
#  50 33 21 50 24 47 42
#  21 33 47 24 42 50 50
#  42 47 21 83 50 15  9
#  50 19 34 42 21 30 24 47
#  24 42 50 21 47 50 33
#  24 47 42 33 50 50 21
#  47 42 21 24 50 33 50
#  24 33 47 42 50 50 21
#  33 50 21 24 42 50 47
#  47 42 24 30 39 34 21 11 19
#  33 42 47 50 50 21 24
#  50 50 42 33 24 21 47
#  33 50 42 68 50 24

set.seed(1)
grp <- find.clusters(gen, max.n.clust = 20, n.pca = 300, choose.n.clust = FALSE, criterion = "goodfit", stat = "BIC")

dapc.xval <- xvalDapc(data_scaled, grp$grp, training.set = 0.9)
    xval.plot.df <- as.data.frame(dapc.xval$DAPC$ind.coord)
    #rownames(xval.plot.df) <- ind_data$Sample_ID_pop
    write_delim(x = cbind(RAD_data2$Sample_ID_pop, xval.plot.df),
                    file = "xval_DAPC_K8.tsv",
                    delim = "\t", col_names = TRUE)
    
     # plot clusters for LDs 1 and 2
    ggplot(data = xval.plot.df, aes(x = LD1, y = LD2, color = grp$grp)) +
    geom_point(cex = 3) +
    theme_classic() +
    labs(x = "Linear discriminant 1", y = "Linear discriminant 2")
    ggsave("LD1-LD2.pdf")

    # plot clusters for next two LDs
    ggplot(data= xval.plot.df, aes(x = LD3, y = LD4, color = grp$grp)) +
    geom_point(cex = 3) +
    theme_classic() +
    labs(x = "Linear discriminant 3", y = "Linear discriminant 4")
    ggsave("LD3-LD4.pdf")