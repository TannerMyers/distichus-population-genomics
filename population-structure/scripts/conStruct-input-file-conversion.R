####################################################################
# This script converts .structure file obtained from Stacks to the #
# format required by the R package 'conStruct'                     #
####################################################################


# Maybe unnecessary see `structure2conStruct` function of 'conStruct' package
# Enter your working directory here
setwd("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/") 

library(tidyverse)


files <- list.files(path="population-structure/conStruct", pattern = "*.structure",full.names=T)
newnames <- paste(sep="",sub('.structure', '',files),"-Processed.structure")

#Loop over all files and make the processed files needed for conStruct
for(i in 1:length(files)){
  
  #Read data file and convert missing data to NA
  str <- read.table(files[1],header=T, sep = "\t", check.names = FALSE) 
    ## Change to 'i' if looping over multiple files
    ## Adding the sep argument prevented the "more columns than column names" error
  # Give first two columns with individual metdata names
  colnames(str)[1:2] <- c("SampleID", "PopmapID")
  # Change missing data to NAs
  str[str == "-9"] <- NA
  
  # Load population identifiers 
  ## This is necessary 
  #write.csv(str[,1:2], file="population-structure/conStruct/structure-file-sample-information.csv", row.names = FALSE)
  PopID <- as.data.frame(read.csv("population-structure/conStruct/structure-file-sample-information.csv", 
                    header=TRUE)[,3])
  colnames(PopID) <- "PopID"
  # Add to `str` object
  str <- cbind(str[,1:2], PopID, str[,3:ncol(str)])
  
  str <- str[ order(str$PopID,str$SampleID),]
  
  #Count number of samples
  SampleID <- as.character(unique(str$SampleID))
  
  #Looping over all loci, create a frequency table of alleles (0,1,2)
  #Jacob Burkhart wrote this loop
  count <- data.frame(SampleID)
  for(loci in 4:dim(str)[2]){   # Changed to '4' from '3' due to addition of another column with population identifiers
    temp <- table(str$SampleID, str[,loci])           
    colnames(temp) <- paste0(colnames(str)[loci], "-", colnames(temp)) 
    temp <- data.frame(unclass(temp)) 
    
    #If there are no alleles, recode the row as -9
    for(j in 1:dim(temp)[1]){
      if(sum(temp[j,]) == 0) {temp[j,] <- NA} 
    }
    #Test if a monomorphic locus slipped through your data processing
    #If so, column bind data to sample ID and any previous datasets
    #If not (as expected), then the column bind will be applied to the 2nd allele
    #Why the 2nd allele?  Because any loci with missing data will result in data being added to the table
    count <- as.matrix(cbind(count,if(length(temp)==1){temp[,1]} else{temp[,2]}))
  }
  
  #Create a vector of the sampling site information for each sample
  pop.vec <- as.vector(str[,3])
  pop.vec <- pop.vec[c(seq(from=1, to=nrow(str), by=2))]
  
  #Make variables to utilize below
  n.pops <- length(unique(pop.vec))
  table.pops <- data.frame(table(pop.vec))
  
  #Make a file of individual sample allele frequencies
  #If you only have one sample per sampling site, then you could stop here
  freq <- matrix(as.numeric(count[,-1])/2,nrow(count),ncol(count)-1)
  f <- matrix(as.numeric(freq),nrow(freq),ncol(freq))
  
  #Empty matrix for sampling site level calculations
  admix.props <- matrix(NA, n.pops,ncol(f))
  
  #Calculate frequency (of 2nd allele) per sampling site
  #The last line tests if there is a sampling site with n=1
  #If so, prints vector because frequency has already been calculated (0, 0.5, or 1)
  #If not, then calculates mean across samples from that site
  for(m in 1:length(table.pops$pop.vec)){
    t<-as.factor(unique(pop.vec))[m]
    admix.props[m,] <- if(table.pops[table.pops$pop.vec == t,2] == 1){f[which(pop.vec==t),]} else{colMeans(f[which(pop.vec==t),],na.rm=T)}
  }
  
  #Export conStruct file and save in working directory
  #write.table(admix.props, file=newnames[1],quote=F,sep="\t",row.names=F,col.names=F)
  write.table(admix.props, file="population-structure/conStruct/populations.structure-Processed.structure",quote=F,sep="\t",row.names=F,col.names=F)
}