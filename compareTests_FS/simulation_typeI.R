# simulation_steadyICA.R  Validate size of the permTest() procedure
# Preliminaries ----
library('steadyICA')
library('MASS')
library('VGAM')
library('parallel')
library('stringr')
source('simulation_functions.R')
set.seed(20210823)

# On Della
#setwd("/scratch/gpfs/eyqian/compareTests")
figPath = 'output/'


# Settings ----
TVec      = c(200, 500, 2000)      # Sample size
df_volVec = c(1, 5, 10, 100, Inf)  # degrees of freedom for shared volatility
nVar      = 3     # Number of variables in DGP
#nSim      = 5000  # Number of simulations
nSim      = 5  # Number of simulations

# Default mixture values for DGP3
delta  = c(0.8 , 1.2 , -1)
kappa  = c(0.06, 0.08, 0.2)
lambda = c(0.52, 0.4 , 0.2)

# Bootstrap settings for tests
#R       = 199
#numboot = 1000


R       = 10
numboot = 10


# Main ----
# Create specification file
Results = list()
jRes    = 1

for(T in TVec){ 
  for(df_vol in df_volVec){
    Results[[jRes]] = list(FUN=simS.DGP1, T=T, nVar=nVar, nSim=nSim, pMat = NA, 
                        args = list(df=5), df_vol = df_vol, 
                        R = R, numboot = numboot)
    Results[[jRes+1]] = list(FUN=simS.DGP2, T=T, nVar=nVar, nSim=nSim, pMat = NA, 
                        args = list(location=0, scale=sqrt(1/2)), df_vol = df_vol,
                        R = R, numboot = numboot)
    Results[[jRes+2]] = list(FUN=simS.DGP3, T=T, nVar=nVar, nSim=nSim, pMat = NA, 
                        args = list(delta = delta, kappa = kappa, lambda = lambda), 
                        df_vol = df_vol, R = R, numboot = numboot)
    jRes = jRes + 3

  }
}

ptm = proc.time()


for(jRes in 1:length(Results))
{
  print(paste(jRes, ' of ', length(Results)))
  # Unpack parameters
  FUN     = Results[[jRes]]$FUN
  T       = Results[[jRes]]$T
  nVar    = Results[[jRes]]$nVar
  args    = Results[[jRes]]$args
  df_vol  = Results[[jRes]]$df_vol  
  R       = Results[[jRes]]$R
  numboot = Results[[jRes]]$numboot
  
  pMat                 = unlist(lapply(1:nSim, function(x) simS(FUN, T, nVar, args, df_vol, R, numboot)))
#  pMat                 = unlist(mclapply(1:nSim, function(x) simS(FUN, T, nVar, args, df_vol, R, numboot)))
  pMat                 = matrix(pMat, nrow = nSim, byrow = TRUE)
  Results[[jRes]]$pMat = pMat
}


time = proc.time() - ptm

# Output results to console ----
# 
# alpha = 0.1
# 
# #T_Results      = unlist(lapply(Results, function(x) x$T))
# df_vol_Results = unlist(lapply(Results, function(x) x$df_vol))
# 
# 
# for(T in TVec){ 
#   for(df_vol in df_volVec){
#     
#     keep = df_vol_Results == df_vol
#     
#     # Setup 10% matrix
#     header              = unlist(lapply(Results, function(x) x$type))
#     mat_size            = matrix(rep(NA, length(Results)*2), nrow = 2)
#     row.names(mat_size) = c( 'DN', 'crossCor')
#     colnames(mat_size)  = header
#     
#     
#     # Add Davis and Ng values
#     temp         = lapply(Results, function(x) x$pMat[,1])
#     pVal         = matrix(unlist(temp), ncol = length(Results))
#     mat_size[1,] = colMeans(pVal <= alpha)
#     
#     # Add cross-correlation values
#     temp   = lapply(Results, function(x) x$pMat[,2])
#     pVal   = matrix(unlist(temp), ncol = length(Results))
#     mat_size[2,] = colMeans(pVal <= alpha)
#   }
# }
# print(round(mat_size,3))

# Save output
save.image(         file = paste0(figPath, 'Results_typeI.RData'))
#write.csv(mat_size, file = paste0(figPath, 'Results_size=10pct.csv'))