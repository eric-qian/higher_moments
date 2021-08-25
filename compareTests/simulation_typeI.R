# simulation_steadyICA.R  Validate size of the permTest() procedure
# Preliminaries ----
library('steadyICA')
library('MASS')
library('parallel')
library('stringr')
source('simulation_functions.R')
set.seed(20210823)

# On Della
setwd("/scratch/gpfs/eyqian/compareTests")
figPath = 'output/'


# Settings
T    = 200
n    = 3
nSim = 5000
pMat = matrix(rep(NA, T*nSim), nrow = nSim)  # Initialize p-value matrix

# Main ----
# Set settings
Results      = list()
Results[[1]] = list(type='normal' , T=T,  n=3,nSim=nSim, pMat = pMat)
Results[[2]] = list(type='t, df=5', T=T, n=3, nSim=nSim, pMat = pMat)
Results[[3]] = list(type='t, df=2', T=T, n=3, nSim=nSim, pMat = pMat)
Results[[4]] = list(type='t, df=1', T=T, n=3, nSim=nSim, pMat = pMat)
Results[[5]] = list(type='ARCH(1)', T=T, n=3, nSim=nSim, pMat = pMat)

ptm = proc.time()
for(jRes in 1:length(Results))
{
  type                 = Results[[jRes]]$type
  nSim                 = Results[[jRes]]$nSim
  pMat                 = unlist(mclapply(1:nSim, function(x) runSim(type, T, n)))
  pMat                 = matrix(pMat, nrow = nSim, byrow = TRUE)
  Results[[jRes]]$pMat = pMat
  print(paste0(type, ': 0.1 nominal rejection Actual = '))
  print(colMeans(pMat < 0.1))
}


time = proc.time() - ptm

# Output results to console ----

# Setup 10% matrix
header              = unlist(lapply(Results, function(x) x$type))
mat_size            = matrix(rep(NA, length(Results)*2), nrow = 2)
row.names(mat_size) = c( 'DN', 'crossCor')
colnames(mat_size)  = header


# Add Davis and Ng values
temp         = lapply(Results, function(x) x$pMat[,1])
pVal         = matrix(unlist(temp), ncol = length(Results))
mat_size[1,] = colMeans(pVal <= .1)

# Add cross-correlation values
temp   = lapply(Results, function(x) x$pMat[,2])
pVal   = matrix(unlist(temp), ncol = length(Results))
mat_size[2,] = colMeans(pVal <= .1)

print(round(mat_size,3))

# Save output
save.image(         file = paste0(figPath, 'Results_typeI.RData'))
write.csv(mat_size, file = paste0(figPath, 'Results_size=10pct.csv'))