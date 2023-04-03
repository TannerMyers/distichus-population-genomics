
library(tidyverse)
library(vcfR)
library(ade4)
library(adegenet)
library(geosphere)
library(LEA)
library(MASS)
library(multcomp)
library(MuMIn)

vcf <- read.vcfR("/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/eems/AnoDist_aln/distichus/distichus-subsp_variants_no-miss_50KbLD.recode.vcf")

RAD_data <- read_table("~/distichus-ddRAD/info/distichus-subsp-Hispaniola-popmap-cleaned.tsv")
colnames(vcf@gt)[2:206] <- RAD_data$Sample_ID_pop

coords <- distinct(RAD_data[, c("Locality", "Longitude", "Latitude")])

data <- vcfR2genind(vcf)
pop(data) <- RAD_data$Locality

pop_data <- genind2genpop(data)

Dgen <- dist.genpop(pop_data, method = 1) # Nei's D

Dgeo <- distm(coords[,2:3])
Dgeo <- as.dist(Dgeo)

ibd <- mantel.randtest(Dgen, Dgeo)
pdf("ibd.pdf")
    plot(ibd)
dev.off()

# taken directly from adegenet tutorial to create 2-D kernel density estimate
pdf("IBD-distichus-subsp.pdf")
    dens <- kde2d(Dgeo, Dgen, n = 300)
    myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
    plot(Dgeo, Dgen, pch = 20, cex = 0.5)
    image(dens, col = transp(myPal(300), 0.7), add = TRUE)
    abline(lm(Dgen ~ Dgeo))
    #title("Isolation by distance plot")
dev.off()

# Now do for north and south paleo island populations
## Note: I am using the same samples I used for hierarchical sNMF analyses, but unfortunately,
## it doesn't appear as though adegenet has an import function for .ped (PLINK format) files ??
## So I used vcftools to create vcfs for these populations, hence the "2" at the end so as not to
## overwrite the log file from creating the plink files
RAD_data_north <- read_table("~/distichus-ddRAD/info/distichus-north-island-popmap-master.tsv")
RAD_data_south <- read_table("~/distichus-ddRAD/info/distichus-south-island-popmap-master.tsv")

n_vcf <- read.vcfR("north-island/north-distichus-ssp2.recode.vcf")
    colnames(n_vcf@gt)[2:176] <- RAD_data_north$Sample_ID_pop
s_vcf <- read.vcfR("south-island/south-distichus-ssp2.recode.vcf")
    colnames(s_vcf@gt)[2:31] <- RAD_data_south$Sample_ID_pop

n_coords <- distinct(RAD_data_north[, c("Locality", "Longitude", "Latitude")])
s_coords <- distinct(RAD_data_south[, c("Locality", "Longitude", "Latitude")])

n_data <- vcfR2genind(n_vcf)
s_data <- vcfR2genind(s_vcf)

pop(n_data) <- RAD_data_north$Locality
pop(s_data) <- RAD_data_south$Locality

n_pop <- genind2genpop(n_data)
s_pop <- genind2genpop(s_data)

n_Dgen<- dist.genpop(n_pop, method = 1)
s_Dgen <- dist.genpop(s_pop, method = 1) # Nei's D

n_Dgeo <- distm(n_coords[, 2:3])
n_Dgeo <- as.dist(n_Dgeo)

s_Dgeo <- distm(s_coords[, 2:3])
s_Dgeo <- as.dist(s_Dgeo)

# north
n_ibd <- mantel.randtest(n_Dgen, n_Dgeo)
pdf("north-island/north_ibd.pdf")
    plot(n_ibd)
dev.off()

# south
s_ibd <- mantel.randtest(s_Dgen, s_Dgeo)
pdf("south-island/south_ibd.pdf")
    plot(s_ibd)
dev.off()

# north
pdf("north-island/IBD-north-distichus-subsp.pdf")
    dens <- kde2d(n_Dgeo, n_Dgen, n = 300)
    myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
    plot(n_Dgeo, n_Dgen, pch = 20, cex = 0.5)
    image(dens, col = transp(myPal(300), 0.7), add = TRUE)
    abline(lm(n_Dgen ~ n_Dgeo))
    #title("Isolation by distance plot")
dev.off()

# south
pdf("south-island/IBD-south-distichus-subsp.pdf")
    dens <- kde2d(s_Dgeo, s_Dgen, n = 300)
    myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
    plot(s_Dgeo, s_Dgen, pch = 20, cex = 0.5)
    image(dens, col = transp(myPal(300), 0.7), add = TRUE)
    abline(lm(s_Dgen ~ s_Dgeo))
    #title("Isolation by distance plot")
dev.off()