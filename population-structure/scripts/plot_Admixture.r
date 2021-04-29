###################################################################
# This script plots Q estimates obtained from Admixture in a bar- #
# plot i.e., the classic structure plot.			              #
###################################################################

# Change to your working directory
working_dir <- getwd()
setwd(working_dir)

# Load library `stringr` for regex 
# install.packages("stringr") # run once
library(stringr)
library(tidyverse)

# Load population map containing individual IDs and lineage information
popmap <- read.table("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/population-structure/admixture/popmap_cleaned_full-lineage-ID.tsv", 
                     sep="\t")
colnames(popmap) <- c("Individual", "Population")

# Create an object containing the different .Q extension files produced by Admixture for the below for loop to 
# iterate over to generate barplots showing ancestry proportions
files <- list.files(#path = "admixture-trial-1/", 
                    pattern = ".Q")

# For each value of K, generate a bar plot to visualize ancestry proportions across individuals in dataset
for (Q in files){
    output <- basename(file.path(Q, fsep=".Q"))
    
    ## Isolate the value of K for the given .Q file
    K <- str_extract(Q, '(\\.[0-9]+\\.)')
    K <- str_extract(K,'[0-9]+') # there is a better way to do this with only 1 line of code, but this works for now

    ## Read in the .Q file output by Admixture as a table
    tbl <- read.table(Q, header=FALSE, fill=TRUE)

    ## Save plots being produced as pdfs
    pdf(file = paste0("Admixture_",Q,".pdf"))
    
    ## Create one dataframe containing Q values and associated individual/population information from the population map
    df <- cbind(popmap, tbl)
    
    ## Re-shape the dataframe such that there are columns for as many values of K and rows for individuals.
    plot_data <- df %>% 
        mutate(id = Individual) %>% #row_number()) %>%
        gather('pop', 'prob', V1:paste0("V",K)) %>% 
        group_by(id) %>%
        mutate(likely_assignment = pop[which.max(prob)], assignment_prob = max(prob)) %>%
        arrange(likely_assignment, desc(assignment_prob)) %>%
        ungroup() %>%
        mutate(id = forcats::fct_inorder(factor(id)))
    
    ## Make barplot of ancestry proportions using ggplot2 & save as object
    admixture_plt <- ggplot(plot_data, aes(id, prob, fill=pop)) + geom_col() + theme_classic()
    
    ## This line of code is necessary to avoid producing empty or corrupted pdf files
    print(admixture_plt)

    dev.off()

}
