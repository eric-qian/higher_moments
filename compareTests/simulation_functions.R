# simulation_functions.R  Functions for Monte Carlo exercise

# Simulate shocks
simS = function(type, T, n)
{
  if(type == 'normal') {  # iid standard normal shocks
    S = mvrnorm(T, rep(0, n), diag(3))
  } else if(type == 't, 5'){  # iid t shocks, df=2
    S = matrix(rt(3*T, 5), ncol=3)
        
  } else if(type == 't, 2'){  # iid t shocks, df=2
    S = matrix(rt(3*T, 2), ncol=3)
    
  } else if(type == 't, 1'){  # iid t shocks, df=1
    S = matrix(rt(3*T, 1), ncol=3)   
    
  } else if (type == 'AR(1), 0.9'){  # Persistent AR(1) shocks
    S     = genAR1(T, n, 0.9)
  } else if (type == 'AR(1), 0.5'){  
    S     = genAR1(T, n, 0.5)
  } else if (type == 'AR(1), 0.1'){
    S     = genAR1(T, n, 0.1)
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
