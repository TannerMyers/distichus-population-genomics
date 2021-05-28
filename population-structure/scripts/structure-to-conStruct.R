
# Load STRUCTURE file as dataframe
str <- read.table("../populations.distichusonly.R0.7.noBahamas.str", header = TRUE, check.names = FALSE)
# Eliminate row of column names
#str <- str[-1,] 

# Recode missing data as NAs
str[str == "-9"] <- NA

#str <- str[ order(str$Population, str$Individual)]

SampleID <- as.character(unique(str$Individual))

# Loop over all loci, creating an allele frequency table (0,1,2)
count <- data.frame(SampleID)
for(loci in 3:dim(str)[2]){   
  temp <- table(str$Individual, str[,loci])           
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
pop.vec <- as.vector(str[,2])
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
write.table(admix.props, "populations.distichusonly.R0.7.noBahamas.processed.str",
            quote=F,
            sep="\t",
            row.names=F,
            col.names=F)