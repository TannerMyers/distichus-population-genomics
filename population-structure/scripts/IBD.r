
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

n.vcf <- read.vcfR("north-island/north-distichus-ssp2.recode.vcf")
    colnames(n.vcf@gt)[2:176] <- RAD_data_north$Sample_ID_pop
s.vcf <- read.vcfR("south-island/south-distichus-ssp2.recode.vcf")
    colnames(s.vcf@gt)[2:31] <- RAD_data_south$Sample_ID_pop

n.coords <- distinct(RAD_data_north[, c("Locality", "Longitude", "Latitude")])
s.coords <- distinct(RAD_data_south[, c("Locality", "Longitude", "Latitude")])

n.data <- vcfR2genind(n.vcf)
s.data <- vcfR2genind(s.vcf)

pop(n.data) <- RAD_data_north$Locality
pop(s.data) <- RAD_data_south$Locality

n.pop <- genind2genpop(n.data)
s.pop <- genind2genpop(s.data)

n.Dgen<- dist.genpop(n.pop, method = 1)
s.Dgen <- dist.genpop(s.pop, method = 1) # Nei's D

n.Dgeo <- distm(n.coords[,2:3])
n.Dgeo <- as.dist(n.Dgeo)

s.Dgeo <- distm(s.coords[,2:3])
s.Dgeo <- as.dist(s.Dgeo)

# north
n.ibd <- mantel.randtest(n.Dgen, n.Dgeo)
pdf("north-island/north_ibd.pdf")
    plot(n.ibd)
dev.off()

# south
s.ibd <- mantel.randtest(s.Dgen, s.Dgeo)
pdf("south-island/south_ibd.pdf")
    plot(s.ibd)
dev.off()

# north
pdf("north-island/IBD-north-distichus-subsp.pdf")
    dens <- kde2d(n.Dgeo, n.Dgen, n = 300)
    myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
    plot(n.Dgeo, n.Dgen, pch = 20, cex = 0.5)
    image(dens, col = transp(myPal(300), 0.7), add = TRUE)
    abline(lm(n.Dgen ~ n.Dgeo))
    #title("Isolation by distance plot")
dev.off() 

# south
pdf("south-island/IBD-south-distichus-subsp.pdf")
    dens <- kde2d(s.Dgeo, s.Dgen, n = 300)
    myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
    plot(s.Dgeo, s.Dgen, pch = 20, cex = 0.5)
    image(dens, col = transp(myPal(300), 0.7), add = TRUE)
    abline(lm(s.Dgen ~ s.Dgeo))
    #title("Isolation by distance plot")
dev.off()