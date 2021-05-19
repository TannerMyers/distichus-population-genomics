#####################################################################
#   This script plots the results of Principal Component Analysis   #
#   performed by plink using the flag `--pca`                       #
#####################################################################

setwd("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/population-structure/pca/plink/") # change to your working directory

###########################################################################

############################# Load libraries #############################
library(tidyverse)

###########################################################################

######################### Load data and clean #############################

# Read in PCA results from plink
pca <- read_table2("populations.r.7.nobrev.pca.eigenvec", col_names = FALSE)
eigenvals <- scan("populations.r.7.nobrev.pca.eigenval")

# Clean up pca object
  ## Remove first individual identifier column. 
  ## Cleaned locality & lineage info will be added shortly for plotting
pca <- pca[,-1] 

  ## Set names for columns. 
  ## The first column contains specimen IDs and the next 20 are principal components
names(pca)[1] <- "Ind"
names(pca)[2:ncol(pca)] <- paste0("PC",1:(ncol(pca)-1))

# Load metadata (lineage and locality) to add to `pca` object
ind_data <- read.csv("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/stacks/info/popmap_cleaned_MacGuigan_no-brev_mapping.csv")[,2:3]

  ## population
pop <- ind_data$Subspecies
  ## locality
loc <- ind_data$REG_Localities

# combine to plot in different colors ## may not be practical for as many as we have
pop_loc <- paste0(pop, "_", loc)

# Remake `pca` dataframe
pca <- as_tibble(data.frame(pca, pop, loc, pop_loc))

###########################################################################

##########################  Plot PCA results ##############################

# Convert to percentage variance explained "PVE"
pve <- data.frame(PC=1:20, pve = eigenvals/sum(eigenvals)*100)

# Create barplot to show the percentage of variance explained by each PC
pve_plot <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
pve_plot + ylab("Percentage variance explained") + theme_light()

# Calculate the cumulative sum of the percent variance explained
cumsum(pve$pve)
#[1]  20.38631  32.74326  41.78656  50.53304  57.10084  62.51664  66.69136  70.29149  73.78208  76.57601
#[11]  79.27787  81.95778  84.46574  86.86103  89.14029  91.39055  93.61343  95.77361  97.91254 100.00000


