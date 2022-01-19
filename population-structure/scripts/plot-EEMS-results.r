setwd('/scratch/tcm0036/distichus-ddRAD/analyses/population-structure/eems/distichus-subsp/data')

library(rEEMSplots)
#library(reemsplots2)

MCMCPATH = c('Hispaniola-distichus-ssp-pruned-EEMS-nDemes-chain1', 'Hispaniola-distichus-ssp-pruned-EEMS-nDemes-chain2', 'Hispaniola-distichus-ssp-pruned-EEMS-nDemes-chain3')
PLOTPATH = 'plots/Hispaniola-distichus-subsp'

eems.plots(mcmcpath = MCMCPATH, plotpath = PLOTPATH, longlat = TRUE, out.png = FALSE)