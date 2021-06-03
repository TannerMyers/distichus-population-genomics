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
pca <- read_table2("populations.cleaned.R0.7.nobrevirostris.eigenvec", col_names = FALSE)
eigenvals <- scan("populations.cleaned.R0.7.nobrevirostris.eigenval")

# Clean up pca object
  ## Remove first individual identifier column. 
  ## Cleaned locality & lineage info will be added shortly for plotting
pca <- pca[,-1] 

  ## Set names for columns. 
  ## The first column contains specimen IDs and the next 20 are principal components
names(pca)[1] <- "Ind"
names(pca)[2:ncol(pca)] <- paste0("PC",1:(ncol(pca)-1))

# Load metadata (lineage and locality) to add to `pca` object
ind_data <- read.csv("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/stacks/info/popmap_cleaned_MacGuigan_no-brev_mapping.csv")

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

# Now plot principal components in a scatter plot

# Make vector with color values 
pop_cols <- c("aurifer" = "royalblue4", "suppar"="royalblue4", "vinosus"="royalblue4", "dom3"="royalblue4", 
              "dom12"="coral4", "dom2"="coral4", 
              "ignigularis"="chocolate1", "igprop"="chocolate1", "ravig"="yellow1", 
              "properus"="grey68", "sejunctus"="grey68", 
              "ravitergum"="yellow1", 
              "dom1"="saddlebrown", "dom4"="saddlebrown", 
              "favillarum"="darkorchid1", 
              "distichus"="palegreen1", "distichoides"="palegreen1", "dapsilis"="palegreen1", "biminiensis"="palegreen1", "ocior"="palegreen1")

# Plot first two principal components
pc_plot <- ggplot(pca, aes(PC1, PC2, col=pop, )) + geom_point(size = 3) #+ geom_text(aes(label=Ind))  
pc_plot <- pc_plot + coord_equal() + theme_light()
pc_plot <- pc_plot + scale_color_manual(values = pop_cols) # Update values here

ggsave(filename = "PC1vPC2inclBahamas.pdf", pc_plot, device = "pdf")

# Now do the same but without Bahamian samples

# Read in PCA results from plink
pca <- read_table2("no-Bahamas/populations.cleaned.R0.7.nobrevirostrisnoBahamas.recoded.eigenvec", 
                   col_names = FALSE)
eigenvals <- scan("no-Bahamas/populations.cleaned.R0.7.nobrevirostrisnoBahamas.recoded.eigenval")

# Clean up pca object
## Remove first individual identifier column. 
## Cleaned locality & lineage info will be added shortly for plotting
pca <- pca[,-1] 

## Set names for columns. 
## The first column contains specimen IDs and the next 20 are principal components
names(pca)[1] <- "Ind"
names(pca)[2:ncol(pca)] <- paste0("PC",1:(ncol(pca)-1))

# Load metadata (lineage and locality) to add to `pca` object
ind_data <- read.csv("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/stacks/info/popmap_cleaned_MacGuigan_no-brev-Bahamas_mapping.csv")

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

# Now plot principal components in a scatter plot

# Make vector with color values 
pop_cols <- c("aurifer" = "royalblue4", "suppar"="royalblue4", "vinosus"="royalblue4", "dom3"="royalblue4", 
              "dom12"="coral4", "dom2"="coral4", 
              "ignigularis"="chocolate1", "igprop"="chocolate1", "ravig"="yellow1", 
              "properus"="grey68", "sejunctus"="grey68", 
              "ravitergum"="yellow1", 
              "dom1"="seagreen4", "dom4"="seagreen4", 
              "favillarum"="darkorchid1"
              )

# Plot first two principal components
pc_plot <- ggplot(pca, aes(PC1, PC2, col=pop, )) + geom_point(size = 3) #+ geom_text(aes(label=Ind))  
pc_plot <- pc_plot + coord_equal() + theme_light()
pc_plot <- pc_plot + scale_color_manual(values = pop_cols) # Update values here

ggsave(filename = "PC1vPC2noBahamas.pdf", pc_plot, device = "pdf")
