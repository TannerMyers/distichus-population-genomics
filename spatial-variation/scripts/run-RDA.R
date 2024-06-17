setwd("/mmfs1/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/gea/rda/LD-pruned/")

library(vcfR)
library(tidyverse)
library(adegenet)
library(vegan)
library(permute)
library(lattice)
library(adespatial)
library(sp)
library(psych)

## Load specimen information
columns <- c('specimen', 'specimen_pop', 'filepath', 'subspecies', 'population', 'complex', 'island', 'locality', 'longitude', 'latitude', 'note')
RAD_data <- read_table("~/distichus-ddRAD/info/distichus-species-group-Hispaniola-only-popmap-cleaned.tsv", col_names = FALSE)
        colnames(RAD_data) <- columns

## Isolate per-sampling locality information
coords <- distinct(RAD_data[,c('locality', 'longitude', 'latitude')])
        dim(coords)

## Load imputed VCF (139,964 SNPs for 191 specimens)
# vcf <- read.vcfR('/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/vcf/distichus-subsp-only/imputation/distichus_geno0.8_ind0.8_maf0.05_dp4_imputed.vcf')
## Load imputed VCF filtered for linkage-disequilibrium (39,427 SNPs for 191 specimens)
vcf <- read.vcfR("/scratch/tcm0036/distichus-ddRAD/alignments/AnoDist/vcf/distichus-subsp-only/imputation/distichus_imp_LD.vcf")
vcf

dat <- vcfR2genind(vcf)
        pop(dat) <- RAD_data$locality

pop_data <- genind2genpop(dat)
pop_data

print(rownames(pop_data@tab) == coords$locality) # populations ordered correctly

## Need to go from counts to frequencies
freq <- makefreq(pop_data)
        ncol(freq)
        nrow(freq)

s <- seq(from = 1, to = ncol(freq), by = 2)
freq_l1 <- freq[, s]

## Perform Hellinger transformation of allele frequencies
freq_h <- decostand(freq_l1, method = "hellinger", MARGIN = 2)

## Load environmental variables
env <- read_table("/scratch/tcm0036/distichus-ddRAD/analyses/AnoDist/selected-environmental-variables.tsv")
env <- distinct(env)
dim(env)
print(env[, c(1:3)] == coords)
env <- env[, -c(1:3)] # drop locality information

print(ls())

## Variable selection, no constraints
# env_only_rda <- rda(freq_h ~ ., data = env, scale = TRUE)
# 	save(env_only_rda, file = "full-environmental-RDA.RData")

# Adjust R^2 for multiple tests
# env_only_r2 <- RsquareAdj(env_only_rda)
# (env_only_r2)

# summary(eigenvals(env_only_rda, model = "constrained"))

# Significance
# sig_env_only <- anova.cca(env_only_rda)
# 	print(sig_env_only)
# 	save(sig_env_only, file = "anova.cca-environment-only.RData")

# Print VIFs
# vif.cca(env_only_rda)

# env_fwd <- forward.sel(freq_h, env, adjR2thresh = env_only_r2, nperm = 9999, alpha = 0.01)
# 	save(env_fwd, file = "environmental-forward-sel.RData")

# variables <- c("clt_min", "elevation_30sec", "rsds_1981.2010_max_V.2.1", "bio15", "hurs_range", "bio3", "bio9", "bio4", "cmi_range")

# fwd_vars <- env[variables]

# sel_env_rda <- rda(freq_h, fwd_vars)
# 	save(sel_env_rda, file = "fwd-selected-environment-RDA.RData")

# sig_sel_env <- anova.cca(sel_env_rda)
# 	save(sig_sel_env, file = "anova-cca-fwdselected-env.RData")

##### Variance partitioning
## Load dbMEMs identified as significant by ``quickMEM" function
load("hellinger-transformed-quickMEM.RData")
dbmems <- quickMEMs$dbMEM_red_model[, -c(2:3)]
print(class(dbmems))
print(dim(dbmems))

## Look at contributions of environment and spatial predictors across fine to broad scales
# vp <- varpart(freq_h, fwd_vars, quickMEMs$dbMEM_red_model[, c(1:2)], quickMEMs$dbMEM_red_model[, c(3:8)])
#	print("Output of 'varpart' performed separating dbMEMs into broad and medium to fine-scale category")
# 	print(vp)

# vp_comb <- varpart(freq_h, fwd_vars, quickMEMs$dbMEM_red_model)
# 	save(vp_comb, file = "var-part-RDA.RData")

#### Variable selection, constraining on space
# env_rda_cond <- rda(freq_h, env, quickMEMs$dbMEM_red_model)
#  	save(env_rda_cond, file = "full-env-dbRDA.RData")
# (env_rda_cond_r2a <- RsquareAdj(env_rda_cond)$r.squared)

# null_m <- rda(freq_h ~ Condition(as.matrix(quickMEMs$dbMEM_red_model)), data = env)

# full_m <- rda(freq_h ~ . + Condition(as.matrix(quickMEMs$dbMEM_red_model)), env)

# ordi <- ordistep(null_m, scope = formula(full_m), direction = "forward")
#  	save(ordi, file = "fwdsel-dbRDA.RData")

## Environmental variables identified by ``ordistep" function when conditioned on space:
# cond_env <- data.frame(env$clt_min, env$elevation_30sec, env$bio4, env$clt_max, env$bio18, env$bio9, env$rsds_1981.2010_max_V.2.1, env$cmi_range, env$bio3)

## CLT_MAX and RSDS_MAX removed because of VIFs over 10
ordi_vars <- c("clt_min", "elevation_30sec", "bio4", "bio18", "bio9", "cmi_range", "bio3")
cond_env <- env[ordi_vars]

# dbRDA <- rda(freq_h, cond_env, quickMEMs$dbMEM_red_model)
#	save(dbRDA, file = "conditional-environment-dbRDA.RData")
# load("conditional-environment-dbRDA.RData")

## Perform dbRDA written in formula interface
dbRDA_scaleF <- rda(freq_h ~ clt_min + elevation_30sec + bio4 + bio18 + bio9 + cmi_range + bio3 + Condition(as.matrix(dbmems)), data = cond_env) #, scale = TRUE)
	# save(dbRDA_2, file = "conditional-environment-dbRDA-2.RData")
	save(dbRDA_scaleF, file = "conditional-environment-dbRDA-scaleF.RData")

# print(summary(dbRDA))
# print(summary(dbRDA_2))
print(summary(dbRDA_scaleF))

## the adjusted R^2 for the dbRDA is
print("The adjusted R^2 for the dbRDA with forward-selected environmental variables is")
# print(RsquareAdj(dbRDA))
#print(RsquareAdj(dbRDA_2))
print(RsquareAdj(dbRDA_scaleF))

## calculate variance inflation factors
# print(vif.cca(dbRDA))
# print(vif.cca(dbRDA_2))
print(vif.cca(dbRDA_scaleF))

# sig_dbRDA <- anova.cca(dbRDA)
#         save(sig_dbRDA, file = "anova-cca-dbRDA.RData")
# sig_dbRDA_2 <- anova.cca(dbRDA_2)
 	# save(sig_dbRDA_2, file = "anova-cca-dbRDA-2.RData")
sig_dbRDA_scaleF <- anova.cca(dbRDA_scaleF)
        save(sig_dbRDA_scaleF, file = "anova-cca-dbRDA-scaleF.RData")


## significance of each RDA axis
# sig_by_axis <- anova.cca(dbRDA_2, by = "axis", parallel = getOption("mc.cores"))
# 	print(sig_by_axis)
# 	save(sig_by_axis, file = "significant-axes-dbRDA-2.RData")
sig_axis_scaleF <- anova.cca(dbRDA_scaleF, by = "axis", parallel = getOption("mc.cores"))
        print(sig_axis_scaleF)
        save(sig_axis_scaleF, file = "significant-axes-dbRDA-scaleF.RData")

