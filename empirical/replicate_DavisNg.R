library(PearsonDS)  # Draw Pearson
library(steadyICA)
library(portes)
library(vars)

# Settings
T    = 400
nSim = 1000

# Matrices from paper
A0    = matrix(c(0.2, 0, 0, 0.3, 0.6, 0, 0.4, 0.3, 0.8), nrow = 3, byrow = TRUE)
B_NLT = matrix(c(1, 3, 0, 1, 1, 0, 0, 0, 1), nrow = 3, byrow = TRUE)
B_LT  = matrix(c(1, 3, 0, 1, 1, 0, 0, 0, 1), nrow = 3, byrow = TRUE)


# Store specifications
Spec = list()
Spec[[1]] = list(Settings = c(BType = 'NLT', uType = 'HL'))
Spec[[2]] = list(Settings = c(BType = 'LT',  uType = 'HL'))
Spec[[3]] = list(Settings = c(BType = 'NLT', uType = 'LL'))
Spec[[4]] = list(Settings = c(BType = 'LT',  uType = 'LL'))


Spec = lapply(Spec, function(x) addParams(x, A0, T, nSim))
Spec = lapply(Spec, function(x) assignB(x, B_NLT, B_LT))


Spec_j = Spec[[1]]


# Functions ---- 

assignB = function(Spec_j, B_NLT, B_LT){
  
  BType = Spec_j$Settings[['BType']]
  uType = Spec_j$Settings[['uType']]
  
  if(BType == 'NLT'){
    B = B_NLT 
  } else if (BType == 'LT'){
      B = B_LT
  }else{
    stop('Enter valid BType setting...')    
  }
  
  Spec_j$B = B  
  return(Spec_j)
}

addParams = function(Spec_j, A, T, nSim){
  Spec_j$A = A 
  Spec_j$T = T
  Spec_j$nSim = nSim
  return(Spec_j)
}


drawLL = function(T)
  {
  moments1 = c(mean=0, variance=1, skewness=2,  kurtosis=20)
  moments2 = c(mean=0, variance=1, skewness=-2, kurtosis=10)
  u1       = rpearson(T, moments=moments1)
  u2       = rpearson(T, moments=moments2)
  u3       = rt(T, df=15)
  u        = cbind(u1, u2, u3)
  return(u)  
} 

drawHL = function(T)
{
  u1       = rStable(T, 1.1, 0)
  u2       = rt(T, df = 5)
  u3       = rt(T, df = 10)
  u        = cbind(u1, u2, u3)
  return(u)  
}

runSim = function(Spec_j, n)
{
  if(Spec_j$Settings['uType'] == 'HL'){
    
    draw_u = match.fun(drawHL)
    
  } else if (Spec_j$Settings['uType'] == 'LL'){
    draw_u = match.fun(drawHL)    
  } else{
    stop('Enter valid shock type')
  }
  
  pValue = data.frame(matrix(data=NA, nrow = Spec_j$nSim, ncol = 11))
  names(pValue) = c('u', 'u2', 'e0', 'etil0', 'etil1', 'etil2', 'etil3',
                     'uhat_etil0', 'uhat_etil1', 'uhat_etil2', 'uhat_etil3')
  

  
  # Simulate
  for(jSim in 1:Spec_j$nSim){
    
    if(jSim %% 10 == 0) print(paste0(jSim, ' of ', nSim, '...'))
    
    # Draw shocks
    u  = draw_u(T)    
    e0 = t(Spec_j$B %*% t(u))
  
    # Simulate data  
    Y     = matrix(NA, nrow = Spec_j$T, ncol = nrow(Spec_j$A))
    Y0    = matrix(0, nrow = nrow(Spec_j$A), ncol = 1)
    Y[1,] = A %*% Y0 + e0[1,]
    for(t in 2:Spec_j$T){
      Y[t,] = A %*% Y[t-1,] + e0[t,]
    }
    
    # Estimate VAR
    res = VAR(Y, 1, type = 'none')
    AHat = Acoef(res)
    eHat = residuals(res)
    
    
    
    pValue$u[jSim] = permTest(u, group = 1:n, R=199, FUN = 'gmultidcov', 
            symmetric = FALSE, alpha=1)
    
    
    
    
  }
  
}


wChol = function(x){
  P = t(chol(var(x)))
  return(t(solve(P) %*% t(x)))
}

wSVD = function(x){
  out = svd(var(x))
  out$v %*% diag(out$d) %*% t(out$u)
  
    P   =t(out$v) %*% diag(out$d)
    
    var(t(P %*% t(x)))
    
  stop('Incomplete')

  prcomp(x)
  
    
  var(t(solve(P) %*% t(x)))
  return(t(solve(P) %*% t(x)))


  }
