###############

# Project: Using ddRADseq data to infer patterns of gene flow across the Anolis distichus
# species complex in Hispaniola

# Authors:
# Tanner Myers, Pietro de Mello, Paul Hime, and Rich Glor

# This script performs cross-validation to compare models under different K 
# values estimated by conStruct

###############

# Change to your working directory
setwd("/gpfs01/home/tcm0036/distichus/population-structure/conStruct/xvalidation")

# Load packages 
## Include the "lib.loc" argument is necessary to tell R on the cluster where to look for the installed package
library(optparse)
library(gtools, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(parallel, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(foreach, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(doParallel, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(Rcpp, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(RcppParallel, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(RcppEigen, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")
library(rstan, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library") 
# the lib.loc argument is necessary to tell R on the cluster where to look for the installed package
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)
#library(devtools)
#install_github("gbradburd/conStruct",build_vignettes=TRUE) # Run once
library(conStruct, lib.loc="/home/tcm0036/miniconda3/envs/R_env/lib/R/library")

# Use processed structure file obtained by running "structure-to-conStruct.R"
freq <- read.table("../populations.distichusonly.R0.7.noBahamas.processed.str")
freq <- as.matrix(freq)

# Load geographic distance matrix generated with topoDistance package
load(file = "../ddRAD_topoDistance_geog_mat.rda") # object's name will be 'tDist.mat'

# Load lat. and long. for sampling localities
latlong <- read.csv("../ddRAD_cleaned_Hispaniolan_distichus_localities.csv", # Update path here
                   header=TRUE)
latlong <- as.matrix(latlong[, c('Longitude', 'Latitude')])
class(latlong)

# 
cl <- makeCluster(8, type="FORK")
registerDoParallel(cl)

my.xvals <- x.validation(train.prop = 0.6, 
                test.prop = 0.4, 
                n.reps = 8,
                K = 1:10, 
                freqs = freq,
                geoDist = tDist.mat,
                coords = latlong,
                prefix = "Hispaniolan_A_distichus",
                data.partitions = NULL, # Return to this post-Evolution meeting as I do have genome position data
                n.iter = 1e4,
                make.figs = FALSE,
                save.files = TRUE,
                parallel = TRUE,
                n.nodes = 8 
                )

stopCluster(cl)

# Read in the results from text files
sp.results <- as.matrix(
                read.table("Hispaniolan_A_distichus_sp_xval_results.txt",
                           header = TRUE,
                           stringsAsFactors = FALSE)
               )
nsp.results <- as.matrix(
                read.table("Hispaniolan_A_distichus_nsp_xval_results.txt",
                           header = TRUE,
                           stringsAsFactors = FALSE)
               )

# Get the 95% confidence intervals for the spatial and nonspatial
#   models over values of K (mean +/- 1.96 the standard error)

sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# then, plot cross-validation results for K=1:10 with 8 replicates

par(mfrow=c(1,2))
pdf(file = "cross_validation_results.pdf")
plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.results,nsp.results),
     main="cross-validation results")
    points(rowMeans(nsp.results),col="green",pch=19)

# finally, visualize results for the spatial model
#   separately with its confidence interval bars
#
# note that you could do the same with the spatial model, 
#   but the confidence intervals don't really show up 
#   because the differences between predictive accuracies
#   across values of K are so large.

plot(rowMeans(sp.results),
     pch=19,col="blue",
     ylab="predictive accuracy",xlab="values of K",
     ylim=range(sp.CIs),
     main="spatial cross-validation results")
segments(x0 = 1:nrow(sp.results),
         y0 = sp.CIs[1,],
         x1 = 1:nrow(sp.results),
         y1 = sp.CIs[2,],
         col = "blue",lwd=2)
dev.off()
