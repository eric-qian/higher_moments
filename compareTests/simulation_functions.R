# simulation_functions.R  Functions for Monte Carlo exercise

# Simulate shocks
simS = function(type, T, n)
{
  if(type == 'normal') {  # iid standard normal shocks
    S = mvrnorm(T, rep(0, n), diag(3))
    
  } else if(str_detect(type, 't, df=')){  # iid t shocks
    df = as.numeric(sapply(strsplit(type, "="), "[", 2))
    S  = matrix(rt(n*T, df), ncol=n)
    
  } else if(str_detect(type, 'AR(1), rho=')){  # AR1
    rho = as.numeric(sapply(strsplit(type, "="), "[", 2))
    S  = genAR1(T, n, rho)
  } else if(type == 'ARCH(1)'){
    S = genARCHq(T, n, 0.5, 0.5)
    
  }else if (str_detect(type, 'SV')) {  # Stochastic volatility 
    # SV, detecting varying tau. Format is "SV, tau=XX."
    
    if(n != 3){
      stop('SV only written for n=3...')
    }
    
    # Extract tau
    tau = as.numeric(sapply(strsplit(type, "="), "[", 2))
    
    # Simulate shocks
    dfs    = c(9, 11, 13)  # t-distribution degrees of freedom
    zeta   = matrix(c(rt(T, dfs[1]), rt(T, dfs[2]), rt(T, dfs[3])), 
                    ncol = 3, byrow = FALSE)
    zeta   = zeta / matrix(rep(sqrt(dfs/(dfs-2)), T), ncol=3, byrow = TRUE)
    sigma  = matrix(rep(exp(-tau^2+tau*rnorm(T)), 3), ncol=3, byrow = FALSE)
    S      = zeta*sigma # shocks
    
  } else{
    stop('Enter valid type')
  }
  return(S)
}

# Run simulation with permutation test
runSim = function(type, T, n, R=199, numboot = 1000){
  S     = simS(type, T, n)  # Generate shocks
  pPerm = permTest(S, group = 1:n, R, FUN = 'gmultidcov', 
                   symmetric = FALSE, alpha=1)
  pCrossCor = crossCorTest(S, numboot)
  return(c(pPerm = pPerm, pCrossCor = pCrossCor))
}

# Make AR1
genAR1 = function(T, n, rho){
  S     = matrix(rep(NA, T*n), nrow = T)
  S0    = rep(0, n)
  S[1,] = rho * S0 + rnorm(n)
  for(t in 2:T){
    S[t,] = rho*S[t-1,] + rnorm(n)
  }
  return(S)
}

# Make mutually independent ARCHq process for process of length T, n variables,
# alpha0, and alphaLagVec = c(alpha1,...,alphaq) for model
#    S_t = sqrt(h_t) * v_t,    v_t iid~ N(0,1)
#    h_t = alpha_0 + alpha_1 S_t-1 ^2 +...+alpha_q S_t-q ^2
genARCHq = function(T, n, alpha0, alphaLagVec){
  q         = length(alphaLagVec)  # ARCH order
  S         = matrix(rep(NA, (T+q)*n), ncol = n)
  S[1:q, ]  = matrix(rep(0, q*n), nrow = q)
  
  for(t in (1+q):nrow(S)){
    
    h_t = rep(alpha0, n)         
    for(l in 1:q){  # Loop through lags
      h_t = h_t + alphaLagVec[l] * S[t-l,]^2
    }
    S[t, ] = sqrt(h_t) * rnorm(n) 
    
  }
  
  S = S[(q+1):nrow(S),]
  return(S)
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
  teststat = crossCorTest.stat(S)
  teststatBoot = sapply(1:numboot, function(x) crossCorTest.boot(S))
  p = mean(teststatBoot > teststat)
  return(p)
  
}

whiten.chol = function(x){
  P = t(chol(var(x)))
  return(t(solve(P) %*% t(x)))
}
