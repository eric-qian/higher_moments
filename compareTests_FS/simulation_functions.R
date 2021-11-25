# simulation_functions.R  Functions for Monte Carlo exercise

# 3 independent, identical t distributions
simS.DGP1 = function(T, nVar, args){
  df          = args[['df']]
  args[['n']] = T*nVar
  S           = do.call(rt, args) / sqrt(df/(df-2))  # Normalize so that sd=1
  S           = matrix(S, nrow=T)
  return(S)
}


# 3 independent, identical Laplace distributions 
simS.DGP2 = function(T, nVar, args){
  args[['n']] = T*nVar
  mu          = args[['location']]
  b           = args[['scale']]
  S           = do.call(rlaplace, args)
  S           = matrix(S, nrow=T)
  return(S)
}


# 3 DLSMN distributions. Chek the parametrization on page 51 of the appendix for 
# Fiorentini and Sentana (2021).
simS.DGP3 = function(T, nVar, args){
  
  delta  = args[['delta']]  # Regulates distance between means
  kappa  = args[['kappa']]  # Ratio of variances
  lambda = args[['lambda']] # Controls mixture
  
  # Validate consistency of arguments
  stopifnot(nVar == length(delta))
  stopifnot(nVar == length(kappa))
  stopifnot(nVar == length(lambda))

  # Store parameters  
  mu_1      = delta*(1-lambda) / sqrt(1+lambda*(1-lambda)*delta^2)
  mu_2      = -lambda / (1-lambda) * mu_1  
  sigma2_1  = 1/((1+ lambda*(1-lambda)*delta^2)*(lambda + (1-lambda)*kappa))
  sigma2_2  = kappa * sigma2_1

  # Output matrix
  S = matrix(data=NA, nrow=T, ncol = nVar)
    
  for(j in 1:nVar){

    # Draw mixture components
    y1   = rnorm(T, mean = mu_1[j], sd = sqrt(sigma2_1[j]))
    y2   = rnorm(T, mean = mu_2[j], sd = sqrt(sigma2_2[j]))
    ind1 = runif(T) < lambda[j]
    
    S[ind1, j]  = y1[ind1]
    S[!ind1, j] = y2[!ind1]
  }
  return(S)

}

# Run simulations. Add shared volatility when df_vol < infinity
simS.vol = function(FUN, T, nVar, args, df_vol=Inf){
  
  # Generate independent shocks
  S = FUN(T, nVar, args)
  
  # Add shared volatility term
  if(df_vol < Inf){
    S = S*rchisq(T, df_vol)/df_vol
  }
  return(S)
}

# Run simulation with permutation test
simS = function(FUN, T, nVar, args, df_vol, R=199, numboot = 1000){
  S     = simS.vol(FUN, T, nVar, args, df_vol)  
  pPerm = permTest(S, group = 1:nVar, R, FUN = 'gmultidcov', 
                   symmetric = FALSE, alpha=1)
  pCrossCor = crossCorTest(S, numboot)
  return(c(pPerm = pPerm, pCrossCor = pCrossCor))
}


crossCorTest.stat = function(S){
  n              = ncol(S)
  corr_shocks_sq = cor(S^2) 
  teststat       = sqrt((sum(as.vector(corr_shocks_sq)^2)-n)/(n^2-n))  # Root mean squared off-diagonal correlation
  return(teststat)  
}

crossCorTest.boot = function(S){
  Sboot = matrix(rep(NA, ncol(S)*nrow(S)), nrow = nrow(S))
  
  # Draw independently across columns
  for(j in 1:3){
    Sboot[,j] = S[sample.int(nrow(S), replace = TRUE), j]
  }
  
  teststat_boot = crossCorTest.stat(Sboot)
  return(teststat_boot)
}


crossCorTest = function(S, numboot){
  teststat     = crossCorTest.stat(S)
  teststatBoot = sapply(1:numboot, function(x) crossCorTest.boot(S))
  p = mean(teststatBoot > teststat)
  return(p)
  
}

whiten.chol = function(x){
  P = t(chol(var(x)))
  return(t(solve(P) %*% t(x)))
}
