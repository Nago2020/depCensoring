library(depCensoring)
library(rvinecopulib)
library(rafalib)
tau = 0.75
Copula = "frank"
Dist.T = "weibull"
Dist.C = "lnorm"
par.T = c(2,1)
par.C = c(1,2)
n=1000
simdata<-TCsim(tau,Copula,Dist.T,Dist.C,par.T,par.C,n)
Y = simdata[[1]]
Delta = simdata[[2]]
fitPar <- ParamCop(Y,Delta,Copula,Dist.T,Dist.C)
