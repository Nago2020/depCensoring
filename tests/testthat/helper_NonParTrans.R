library(depCensoring)
library(MASS)
library(survival)
library(pbivnorm)
library(nleqslv)

# data example

rho = 0.8
n = 300
beta = c(0.5, 1); eta = c(1,1.5)
sigma = matrix(c(1,rho,rho,1),ncol=2)
err = MASS::mvrnorm(n, mu = c(0,0) , Sigma=sigma)
err1 = err[,1]; err2 = err[,2]
x1 = rbinom(n,1,0.5); x2 = runif(n,-1,1)
X = matrix(c(x1,x2),ncol=2,nrow=n); W = X   # data matrix
T1 = X%*%beta+err1
C =  W%*%eta+err2
T1 = exp(T1); C = exp(C)
 A = runif(n,0,8); Y = pmin(T1,C,A)
d1 = as.numeric(Y==T1)
d2 = as.numeric(Y==C)
resData = data.frame("Z" = Y,"d1" = d1, "d2" = d2)   # should be data frame
colnames(X) = c("X1", "X2")
colnames(W) = c("W1","W2")

#  Bootstrap is False by default
output = NonParTrans(resData = resData, X = X, W = W, n.iter = 3)
l = length(output$parameterEstimates)
rhohat <- output$parameterEstimates[l]
