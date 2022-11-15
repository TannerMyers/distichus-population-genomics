setwd("/scratch/tcm0036/distichus-ddRAD/analyses/environmental-association/")

rm(list = ls())

# Load R packages
library(sp)
library(raster)
library(rgeos)
library(rgdal)
library(RStoolbox)
library(adegenet)
library(vcfR)
library(hierfstat)
library(StAMPP)
library(BEDASSLE)
library(ecodist) # MRM function may be able to replace MMRR
library(gdm)
library(vegan)
library(tidyverse)

# ddRADseq individual information
RAD_data <- read_table("/home/tcm0036/distichus-ddRAD/info/distichus-popmap-cleaned-master.tsv", col_names = TRUE)
snmfK9 <- read_table("/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/lea/sNMF_K9_ancestry_coefficients_cluster-specimen-cleaned.tsv")
RAD_data2 <- RAD_data %>% filter(Sample != "4481" & Sample != "10085" & Sample != "10086") %>%
    filter(!(Taxon %in% c("websteri", "marron", "caudalis", "altavelensis"))) %>%
    filter(Island %in% "Hispaniola")  %>%
    filter(!grepl("_rep", Sample)) # Temporary, eventually I need to deal with the technical replicates in the dataset

RAD_data2 <- RAD_data2[RAD_data2$Sample_ID_pop %in% snmfK9$V1[snmfK9$cluster != 6 & snmfK9$cluster != 8], ]

# Loop over individual chromosome vcf files
for (i in list.files("vcfs-by-chrom/", pattern = ".vcf")){
    vcfR <- read.vcfR(paste0("vcfs-by-chrom/", i))

    # convert vcfR object to genind object
    dat <- vcfR2genind(vcfR)
    dat2 <- dat[indNames(dat) %in% RAD_data2$Sample_ID_pop]
    pop(dat2) <- RAD_data2$Locality

    data <- genind2hierfstat(dat = dat2, pop = dat2@pop)
    data <- as.data.frame(cbind(dat2@pop, data[, c(2:ncol(data))]))

    gen_diff <- pairwise.WCfst(dat = data, diploid = TRUE)
    save(gen_diff, file = paste0("vcfs-by-chrom/", str_replace(i, ".vcf", ""), "-Fst.RData"))
}