################################################################################
# This script plots Q estimates obtained from Admixture using the R package    #
# `pophelper` v. 2.3.1                                                         #
#################################################################################

# change to your working directory
setwd("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/")

###########################################################################

############################# Load Libraries #############################
library(pophelper) # run `remotes::install_github('royfrancis/pophelper')` first

# check if version is 2.3.1
packageDescription("pophelper", fields="Version")

###########################################################################

################### Load .Q files from Admixture ##########################

# Load population map containing individual IDs and lineage information
popmap <- read.csv("stacks/info/popmap_cleaned_MacGuigan_no-brev_mapping.csv", 
                   stringsAsFactors = FALSE)[,1:3]
colnames(popmap) <- c("Individual", "type", "loc") # Add column names 

sapply(popmap, is.character) # check if character data type

# List files in working directory that have the .Q extension
afiles <- list.files(path="population-structure/admixture/admixture-trial2_no-brev", 
                    pattern=".Q", 
                    full.names= TRUE)

# Convert q-matrices to qlist
alist <- readQ(files = afiles)
lapply(alist, attributes)

# Since all are of equal length, add specimen labels for each individual
if(length(unique(sapply(alist,nrow)))==1) alist <- lapply(alist,"rownames<-",popmap$Individual)

#tabulateQ(qlist = alist, writetable = TRUE, exportpath =getwd()) # doesn't provide new information

###########################################################################

#########################  Plot Admixture runs ############################

output_dir="/population-structure/admixture/admixture-trial2_no-brev/plots/pophelper"

pops <- popmap[,2,drop=FALSE]

plotQ(alist[c(2)], 
      sortind = "all", # Order by clusters
      showindlab = TRUE,
      useindlab = TRUE,
      basesize = 10,
      #grplab = pops, indlabwithgrplab = TRUE,
      #grplabsize = 4, grplabface = "italic",
      #linesize = 0.8, pointsize = 3,
      returnplot = TRUE,
      #exportplot = FALSE) 
      imgtype = "pdf", # set export file as pdf
      exportpath = output_dir) # exports pdfs to this directory


