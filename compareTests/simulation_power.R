# Preliminaries ---
library('steadyICA')
library('MASS')
library('parallel')
library('stringr')
source('simulation_functions.R')
set.seed(2021)

# On Della
setwd("/scratch/gpfs/eyqian/compareTests")
figPath = 'output/'

# Settings
T    = 200
n    = 3
nSim = 5000
pMat = NA
tauVec = seq(0, 1, by=.025)

# Main ----
# Set settings
Results      = list()
for(jRes in 1:length(tauVec)){
  Results[[jRes]] = list(type=paste0('SV, tau=', tauVec[jRes]), T=T, n=3, 
                         nSim=nSim, pMat = pMat)
}


ptm = proc.time()
for(jRes in 1:length(Results))
{
  print(paste('Running tau =', tauVec[jRes], '...'))
  type                 = Results[[jRes]]$type
  nSim                 = Results[[jRes]]$nSim
  pMat                 = unlist(mclapply(1:nSim, function(x) runSim(type, T, n)))
  pMat                 = matrix(pMat, nrow = nSim, byrow = TRUE)
  Results[[jRes]]$pMat = pMat  
}
time = proc.time() - ptm



# Make 10% table ----
rejRate       = data.frame(matrix(rep(NA, length(Results) * 2), 
                                  nrow=length(Results)))
names(rejRate) = c('perm', 'crossCor')

for(i in 1:2){
  temp        = lapply(Results, function(x) x$pMat[,i])
  pVal        = matrix(unlist(temp), ncol = length(Results))
  rejRate[,i] = colMeans(pVal < 0.1)
}

# Export p-values ----
# Write p values DN
temp           = lapply(Results, function(x) x$pMat[,1])
pVal           = matrix(unlist(temp), ncol = length(Results))
header         = unlist(lapply(Results, function(x) x$type))
pVal_df        = data.frame(pVal)
names(pVal_df) = sapply(Results, function(x) x$type)
write.csv(pVal_df, paste0(figPath, 'power_pVal_DN.csv'))

# Write p-values crossCor
temp           = lapply(Results, function(x) x$pMat[,2])
pVal           = matrix(unlist(temp), ncol = length(Results))
header         = unlist(lapply(Results, function(x) x$type))
pVal_df        = data.frame(pVal)
names(pVal_df) = sapply(Results, function(x) x$type)
write.csv(pVal_df, paste0(figPath, 'power_pVal_crossCor.csv'))

# Save files
save.image(file=paste0(figPath, 'Results_power.RData'))
