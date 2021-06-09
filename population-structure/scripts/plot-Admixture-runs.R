
setwd("/gpfs01/home/tcm0036/distichus/population-structure/admixture") # change to your working directory

# Load libraries
library(stringr)
library(dplyr)
library(conStruct)


# Assign variables
## Load file with individual metadata. This will be used to sort Q files and cluster 
## individuals by population
popmap <- read.table("popmap_cleaned_MacGuigan_no-brev_mapping.csv", sep=",",header=TRUE)

## "files" will be looped over, creating a structure plot for each value of Q
files <- list.files(pattern=".Q", recursive=TRUE) # "recursive" argument searches within subdirectories
files <- files[-c(1:15)] # exclude files from preliminary Admixture run with 

# For each value of K, generate a bar plot to visualize ancestry proportions across individuals in dataset
for (Q in files){
  
  # Load name of Q file
  output <- basename(file.path(Q, fsep=".Q"))
  
  # Isolate the Admixture run
  run <- str_extract(output, 'run[0-9]+')
  
  ## Isolate the value of K for the given .Q file
  K <- str_extract(Q, '(\\.[0-9]+\\.Q)') # extract the value of K from the filename
  K <- str_extract(K,'[0-9]+') # strip away everything except the number
  # there is probably a better way to do this with only 1 line of code, but this works for now
  
  ## Read in the .Q file output by Admixture as a table
  df <- read.table(Q, header=FALSE, fill=TRUE)
 
  ## Combine population map
  df2 <- cbind(popmap, df)
  
  # Sort rows by population identity
  df3 <- df2 %>% arrange(K9Cluster)
  
  ## Get dataframe of sorted admixture proportions
  admix.props <- as.matrix(df3[,9:ncol(df3)])
  
  pdf(paste0("plots/",run,".",K,"barplot.pdf"))
  ## Make STRUCTURE barplots using the "make.structure.plot" from conStruct
  make.structure.plot(admix.props, sample.names = df3$PopID, mar = c(5,4,1.5,1.5))
  dev.off()
}
