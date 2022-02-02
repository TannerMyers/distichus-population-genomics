#####################################################################
#   This script plots the results of Principal Component Analysis   #
#   performed by plink using the flag `--pca`                       #
#####################################################################

working_dir <- getwd()
setwd(working_dir) # change to your working directory

###########################################################################

############################# Load libraries #############################
library(tidyverse)

###########################################################################

######################### Load data and clean #############################

# Read in PCA results from plink
#pca <- read_table("population-structure/pca/ref-alignment/variants_minQ20minDP10maxDPnanmac3geno95ind25_FIL-4.eigenvec", col_names = FALSE)
#eigenvals <- read_table("population-structure/pca/ref-alignment/variants_minQ20minDP10maxDPnanmac3geno95ind25_FIL-4.eigenval", col_names = FALSE)

pca <- read_table("population-structure/pca/ref-alignment/linkage-pruned/distichus-only-LD-pruned.eigenvec", col_names = FALSE)

## The second column, which would normally contain locality information instead repeats individual identifiers, 
## except in the case of technical replicate samples. The parser did something odd there, recognizing the '_' between 
## sample ID and 'rep' as a break so the entries with the second column just read as the rest of the file path beginning with 'rep'

#pca <- read_table2("populations.cleaned.R0.7.nobrevirostris.eigenvec", col_names = FALSE) # For plink output of Stacks' assembled loci
#eigenvals <- scan("populations.cleaned.R0.7.nobrevirostris.eigenval")

# Clean up 'pca' object

## Cleaning up the superfluous filepath information used in the popmap supplied to the variant calling script, leaving just sample IDs
pca$X1 <- str_replace(pca$X1, "/scratch/tcm0036/distichus-ddRAD/alignments/results/bam/", "") %>% 
          str_replace(".aligned.sorted.bam","")
## Because the parser split up the filepaths for the techincal replicate samples, do this to reflect which ones are technical replicates
pca$X1 <- with(pca, ifelse(X2=="rep.aligned.sorted.bam", paste0(X1,"_rep"), X1))

## Remove first individual identifier column. 
## Cleaned locality & lineage info will be added shortly for plotting
#pca <- pca[,-1] 
pca <- pca[,-2]

  ## Set names for columns. 
  ## The first column contains specimen IDs and the next 20 are principal components
names(pca)[1] <- "Ind"
names(pca)[2:ncol(pca)] <- paste0("PC",1:(ncol(pca)-1))

# Load metadata (lineage and locality) to add to `pca` object
#ind_data <- read.csv("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/stacks/info/popmap_cleaned_MacGuigan_no-brev_mapping.csv")
ind_data <- read_table('~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/info/distichus-popmap-master.tsv', col_names=T)
## population
#pop <- ind_data$Subspecies
  ## locality
#loc <- ind_data$REG_Localities
## repeat steps to clean sample IDs
dropped <- !ind_data$Sample%in%pca$Ind
ind_data[dropped,] 
ind_data <- ind_data[!dropped,]
pop <- ind_data$PopID
taxon <- ind_data$Taxon
loc <- ind_data$Locality

# combine to plot in different colors ## may not be practical for as many as we have
pop_loc <- paste0(pop, "_", loc)

# Remake `pca` dataframe
pca2 <- as_tibble(data.frame(pca$Ind, pop, loc, pop_loc, pca[,2:ncol(pca)]))

###########################################################################

##########################  Plot PCA results ##############################

# Convert to percentage variance explained "PVE"
pve <- data.frame(PC=1:20, pve = eigenvals/sum(eigenvals)*100)
names(pve)[2] <- "pve"

# Create barplot to show the percentage of variance explained by each PC
pve_plot <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
pve_plot + ylab("Percentage variance explained") + theme_light()

# Calculate the cumulative sum of the percent variance explained
cumsum(pve$pve)

# Now plot principal components in a scatter plot

# Make vector with color values 
pop_cols <- c("distichus"="palegreen1", "distichoides"="palegreen1", "dapsilis"="palegreen1", "biminiensis"="palegreen1", "ocior"="palegreen1",
              "dom1"="tan4", "dom2"="sienna", "dom4"="orange4", 
              "Cordillera_ignigularis"="darkorange3", "west_ignigularis"="orange", "east_ignigularis"="darkorange",
              "properus"="grey68", "sejunctus"="grey68", "igprop"="darkkhaki", "ravig"="khaki1", "ravitergum"="gold1", "altavelensis"="gold1",
              "orange_favillarum"="slateblue1", "white_favillarum"="slateblue4",
              "aurifer" = "royalblue4", "suppar"="royalblue4", "vinosus"="royalblue4", "dom3"="royalblue4")

tax_cols <- c("aurifer" = "royalblue4", "suppar"="royalblue4", "vinosus"="royalblue4", "dom3"="royalblue4", 
              "dom12"="coral4", "dom2"="firebrick", 
              "ignigularis"="chocolate1", "igprop"="chocolate1", "ravig"="yellow1", 
              "properus"="grey68", "sejunctus"="grey68", 
              "ravitergum"="yellow1", 
              "dom1"="saddlebrown", "dom4"="saddlebrown", 
              "favillarum"="darkorchid1", 
              "distichus"="palegreen1", "distichoides"="palegreen1", "dapsilis"="palegreen1", "biminiensis"="palegreen1", "ocior"="palegreen1",
              "altavelensis"="gray0", "brevirostris"="gray0", "caudalis"="gray0", "marron"="gray0", "websteri"="gray0")


# Plot first two principal components
pc_plot <- ggplot(pca2, aes(PC1, PC2, color=pop)) + geom_point(size = 3) #+ geom_text(aes(label=pca.Ind))#, aes(shape = pop_shps)) 
pc_plot <- pc_plot + coord_equal() + theme_light()
pc_plot <- pc_plot + scale_color_manual(values = pop_cols) 
pc_plot

pc_plot_23 <- ggplot(pca2, aes(PC2, PC3, col=pop, )) + geom_point(size = 3)
pc_plot_23 <- pc_plot_23 + coord_equal() + theme_light()
pc_plot_23 <- pc_plot_23 + scale_color_manual(values = pop_cols) 
pc_plot_23


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
pc_plot
ggsave(filename = "PC1vPC2noBahamas.pdf", pc_plot, device = "pdf")
