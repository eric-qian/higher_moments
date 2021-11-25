library('steadyICA')
library('MASS')
library('VGAM')
library('parallel')
library('stringr')
source('simulation_functions.R')
set.seed(20210823)


T    = 100000
nVar = 3

# Check first and second moments of DGPs. Should be mean zero and variance 1. ----
S1 = simS.DGP1(T, nVar, list(df=5))
colMeans(S1)
cov(S1)

S2 = simS.DGP2(T, nVar, list(location = 0, scale=sqrt(1/2)))
colMeans(S2)
cov(S2)

delta  = c(0.8, 1.2, -1)
kappa  = c(0.06, 0.08, 0.2)
lambda = c(0.52, 0.4, 0.2)

S3 = simS.DGP3(T, nVar, args = list(delta = delta, kappa = kappa, lambda = lambda))
colMeans(S3)
cov(S3)

df = 100000
var(rchisq(T, df)/df)