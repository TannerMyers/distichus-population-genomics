setwd("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/population-structure/admixture/R0.7-nobrevirostris-admixture-results/")

# Load library `stringr` for regex 
# install.packages("stringr") # run once
library(stringr)
library(tidyverse)
library(conStruct)
library(maps)

# Load population map containing individual IDs and lineage information
popmap <- read.table(#"~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/stacks/info/popmap_cleaned_nobrev_Admixture.tsv", 
                    "~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/stacks/info/popmap_cleaned_MacGuigan_no-brev_mapping.csv", 
                    sep=",", header=TRUE)
#colnames(popmap) <- c("Individual", "Population", "Subspecies")

# Load .Q file output from Admixture
test <- read.table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/population-structure/admixture/R0.7-nobrevirostris-admixture-results/populations.cleaned.R0.7.nobrevirostris.recoded.9.Q")

# 
test2 <- cbind(popmap, test)

# Sort rows by population identity
test3 <- test2 %>% arrange(PopID)

test4 <- test3[test3$Island=="Hispaniola",] # to map only Hispaniolan taxa

admix.props <- as.matrix(test3[,8:ncol(test3)])



# Use the "make.structure.plot" function from conStruct to plot admixture proportions
pdf("K9.Admixprops.pdf")
make.structure.plot(admix.props, sample.names = test3$PopID, mar = c(5,4,1.5,1.5))
dev.off()

pdf("K9_pieplot.pdf")
# Need to include dataframe with coords in long, lat (x, y) order for sampling localities
#make.admix.pie.plot(admix.proportions = admix.props, coords = test3[,c(6,5)], radii=1)

maps::map(xlim = range(test4$Longitude) 
          +c(-1,1), 
           ylim = range(test4$Latitude)
          +c(-1,1), 
           col="gray")
make.admix.pie.plot(admix.proportions = admix.props,
                    coords = test3[,c(6,5)],
                    add = TRUE, radii = 1)

dev.off()