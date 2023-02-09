
working_dir <- getwd()
setwd(working_dir)

library(tidyverse)

ind_data <- read_table('~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/info/distichus-popmap-master.tsv', col_names=T)

pca <- read_table(list.files(pattern = "eigenvec"), col_names = FALSE)
eigenvals <- read_table(list.files(pattern = "eigenval"), col_names = FALSE)

# Cleaning up the superfluous filepath information used in the popmap supplied to the variant calling script, leaving just sample IDs
pca$X1 <- str_replace(pca$X1, "/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/bam/", "") %>%
          str_replace(".aligned.sorted.bam", "")
# Because the parser split up the filepaths for the techincal replicate samples, do this to reflect which ones are technical replicates
pca$X1 <- with(pca, ifelse(X2 == "rep.aligned.sorted.bam", paste0(X1, "_rep"), X1))

## Remove first individual identifier column.
## Cleaned locality & lineage info will be added shortly for plotting
pca <- pca[, -2]

names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

dropped <- !ind_data$Sample %in% pca$ind
ind_data[dropped,] 
ind_data <- ind_data[!dropped,]
pop <- ind_data$PopID
taxon <- ind_data$Taxon
loc <- ind_data$Locality
# combine to plot in different colors ## may not be practical for as many as we have
pop_loc <- paste0(pop, "_", loc)

# Remake `pca` dataframe
pca2 <- as_tibble(data.frame(pca$ind, pop, taxon, loc, pop_loc, pca[, 2:ncol(pca)]))

tax_cols <- c("aurifer" = "royalblue4", "suppar"="royalblue4", "vinosus"="royalblue4", "dom3"="royalblue4",
              "dom12" = "coral4", "dom2" = "firebrick",
              "ignigularis" = "chocolate1", "igprop" = "chocolate1", "ravig" = "yellow1",
              "properus" = "grey68", "sejunctus" = "grey68",
              "ravitergum" = "yellow1", "altavelensis" = "yellow1",
              "dom1" = "saddlebrown", "dom4" = "saddlebrown",
              "favillarum" = "darkorchid1",
              "distichus" = "palegreen1", "distichoides" = "palegreen1", "dapsilis" = "palegreen1", "biminiensis" = "palegreen1", "ocior" = "palegreen1",
              "brevirostris" = "gray0", "caudalis" = "gray0", "marron" = "gray0", "websteri" = "gray0")

# Convert to percentage variance explained "PVE"
pve <- data.frame(PC = 1:20, pve = eigenvals/sum(eigenvals)*100)
    names(pve)[2] <- "pve"

# Create barplot to show the percentage of variance explained by each PC
pve_plot <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
pve_plot + ylab("Percentage variance explained") + theme_light()

# Calculate the cumulative sum of the percent variance explained
cumsum(pve$pve)

# Plot first two principal components
pc_plot <- ggplot(pca2, aes(PC1, PC2, color = taxon)) + geom_point(size = 3) #+ geom_text(aes(label=pca.ind))#, aes(shape = pop_shps)) 
pc_plot <- pc_plot + coord_equal() + theme_light()
pc_plot <- pc_plot + scale_color_manual(values = tax_cols)
pc_plot

pc_plot_23 <- ggplot(pca2, aes(PC2, PC3, color = taxon)) + geom_point(size = 3)
pc_plot_23 <- pc_plot_23 + coord_equal() + theme_light()
pc_plot_23 <- pc_plot_23 + scale_color_manual(values = tax_cols)
pc_plot_23

#########################################################################################################
# ADMIXTURE

files <- list.files(pattern = "[0-9].Q", recursive = TRUE)

for (Q in files){
    # Load name of Q file
    input <- basename(file.path(Q, fsep = ".Q"))

    ## Isolate the value of K for the given .Q file
    K <- str_extract(Q, '(\\.[0-9]+\\.Q)') %>% str_extract('[0-9]+')

    ## Read in the .Q file output by Admixture as a table
    df <- read.table(Q, header = FALSE, fill = TRUE)

    ## Combine population map
    df2 <- cbind(ind_data, df)

    # Sort rows by population identity
    df3 <- df2 %>% arrange(Taxon) # arrange(Taxon)

    ## Get dataframe of sorted admixture proportions
    #admix.props <- as.matrix(df3[,9:ncol(df3)])
    admix.props <- as.matrix(df3[, 13:ncol(df3)])

    pdf(paste0("barplot-", K, ".pdf"), height = 7, width = 35)
        ## Make STRUCTURE barplots using the "make.structure.plot" from conStruct
        make.structure.plot(admix.props, sample.names = df3$Sample_ID_pop, mar = c(5,4,1.5,1.5))
    dev.off()

    #maps::map(xlim = range(popmap$Longitude)
    #        +c(-1,1), 
    #        ylim = range(popmap$Latitude)
    #        +c(-1,1), 
    #        col="gray")
    #pdf("pieplot-",K,".pdf")
    #make.admix.pie.plot(admix.proportions = admix.props,
    #                  coords = df3[,c(7,8)],
    #                  add = TRUE, radii = 1)
    #dev.off()
}
