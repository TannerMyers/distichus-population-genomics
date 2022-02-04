library(foreach)
library(doMC)
library(pegas)

# Load SNPs2CF function from Melissa Olave's GitHub: https://github.com/melisaolave/SNPs2CF
source("/home/tcm0036/distichus-ddRAD/scripts/SNPs2CF_functions.R") 

# Provide matrix of SNP data (e.g., a phylip file) and file with
# individual - lineage assignmments as command line arguments 
# For example, "Rscript variants.phy"
args <- commandArgs(trailingOnly = TRUE)   
seqMatrix <- args[1]
IMap <- args[2]

SNPs2CF(wd = getwd(),
        seqMatrix = seqMatrix, 
        ImapName = IMap,
        n.quartets = 3,
        max.SNPs = NULL,
        between.sp.only = TRUE,
        bootstrap = TRUE,
        save.progress=TRUE,
        cores = 4)