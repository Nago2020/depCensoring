% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NonparametricTransformation.R
\name{NonParTrans}
\alias{NonParTrans}
\title{Fit a semiparametric transformation model for dependent censoring}
\usage{
NonParTrans(
  resData,
  X,
  W,
  start = NULL,
  n.iter = 15,
  bootstrap = FALSE,
  n.boot = 50,
  eps = 0.001
)
}
\arguments{
\item{resData}{Data matrix with three columns;  Z = the observed survival time, d1 = the censoring indicator of T
and  d2 =  the censoring indicator of C.}

\item{X}{Data matrix with covariates related to T}

\item{W}{Data matrix with covariates related to C}

\item{start}{Initial values for the finite dimensional parameters. If \code{start} is NULL, the initial values will be obtained
by fitting an Accelerated failure time models.}

\item{n.iter}{Number of iterations; the default is \code{n.iter = 20}. The larger the number of iterations, the longer the computational time.}

\item{bootstrap}{A boolean indicating whether to compute bootstrap standard errors for making inferences.}

\item{n.boot}{Number of bootstrap samples to use in the estimation of bootstrap standard errors if \code{bootstrap = TRUE}. The default is n.boot = 50. But, higher
values  of \code{n.boot} are recommended for obtaining good estimates of bootstrap standard errors.}

\item{eps}{Convergence error. This is set by the user in such away that the desired convergence is met; the default is \code{eps = 1e-3}.}
}
\value{
This function returns a fit of a semiparametric transformation model; parameter estimates, estimate of the non-parametric transformation function, bootstrap standard
errors for finite-dimensional parameters, the nonparametric cumulative hazard function, etc.
}
\description{
This function allows to estimate the dependency parameter along all other model parameters. First, estimates a non-parametric transformation function, and
then at the second stage it estimates other model parameters assuming that the non-parametric function is known. The details for
implementing the dependent censoring methodology can be found in Deresa and Van Keilegom (2021).
}
\examples{
\donttest{
# Toy data example to illustrate implementation
n = 300
beta = c(0.5, 1); eta = c(1,1.5); rho = 0.70
sigma = matrix(c(1,rho,rho,1),ncol=2)
err = MASS::mvrnorm(n, mu = c(0,0) , Sigma=sigma)
err1 = err[,1]; err2 = err[,2]
x1 = rbinom(n,1,0.5); x2 = runif(n,-1,1)
X = matrix(c(x1,x2),ncol=2,nrow=n); W = X   # data matrix
T1 = X\%*\%beta+err1
C =  W\%*\%eta+err2
T1 = exp(T1); C = exp(C)
A = runif(n,0,8); Y = pmin(T1,C,A)
d1 = as.numeric(Y==T1)
d2 = as.numeric(Y==C)
resData = data.frame("Z" = Y,"d1" = d1, "d2" = d2)   # should be data frame
colnames(X) = c("X1", "X2")
colnames(W) = c("W1","W2")

#  Bootstrap is false by default
output = NonParTrans(resData = resData, X = X, W = W, n.iter = 2)
output$parameterEstimates

}
}
\references{
Deresa, N. and Van Keilegom, I. (2021). On semiparametric modelling, estimation and inference for survival data subject to dependent censoring, Biometrika, 108, 965–979.
}
