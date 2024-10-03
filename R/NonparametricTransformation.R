#' @title Estimate a nonparametric transformation function
#'
#' @description This function estimates the nonparametric transformation  function H when the survival time and censoring time
#' are dependent given covariates.  The estimating equation of H was derived based on the martingale ideas. More details about
#' the derivation of a nonparmaetric estimator of H and its estimation algorithm can be found in Deresa and Van Keilegom (2021).
#'
#' @references  Deresa, N. and Van Keilegom, I. (2021). On semiparametric modelling, estimation and inference for survival data subject to dependent censoring, Biometrika, 108, 965–979.
#'
#' @param theta Vector of parameters in the semiparametric transformation model.
#' @param resData Data matrix with three columns;  Z = the observed survival time, d1 = the censoring indicator of T
#'   and  d2 =  the censoring indicator of C.
#' @param X Data matrix with covariates related to T.
#' @param W Data matrix with covariates related to C.
#' @importFrom stats pnorm  dnorm qnorm sd
#' @import pbivnorm
#' @import nleqslv
#'
#' @return Returns the estimated transformation function H for a fixed value of parameters theta.
#'
#' @export


SolveH <- function(theta,resData, X, W){                                 # Z,nu,X,W,beta,eta,rho
  k = ncol(X)
  l = ncol(W)
  beta = theta[1:k]
  eta = theta[(k+1):(l+k)]
  rho = theta[(l+k+1)]
  Z = resData$Z
  d1 = resData$d1
  d2 = resData$d2
  n = length(Z)
  nu = d1+d2
  T1 = c(-Inf,Z[nu==1])
  T1 = unique(T1)                                                     # distinct observed values of Y
  m = length(T1);
  H = rep(0,m);
  H[1] = -Inf;                                                        # for the first observed time, H goes to negative infinity by assumption
  # H[2] = uniroot(SolveHt1,interval = c(-7,7),Z=Z,nu=nu,t=T1[2],X=X,W=W,beta=beta,eta=eta,rho=rho)$root
  H[2] = nleqslv(0,SolveHt1,jac = NULL,Z=Z,nu=nu,t=T1[2],X=X,W=W,theta)$x;

  for (k in 3:m)
  {
    Sum = 0
    j = k-1
    for (i in 1:n)
    {
      upx = H[j]-X[i,]%*%beta
      upy = H[j]-W[i,]%*%eta
      bvprob = Bvprob(upx,upy,rho);
      der = (-dnorm(upy)*(1-pnorm((upx-rho*upy)/(1-rho^2)^0.5))-dnorm(upx)*(1-pnorm((upy-rho*upx)/(1-rho^2)^0.5)))
      if (bvprob!=0)
        Sum = Sum +(Z[i]>=T1[k])/bvprob*der
      else
        Sum = Sum
    }
    if (Sum!=0)
    {
      H[k] = H[j]-sum((Z==T1[k])*(nu==1))/Sum
    }else {
      H[k] = H[j]
    }
  }
  return(cbind(T1,H))
}

#' @title Estimating equation for Ht1
#'
#' @description This function obtains an estimating equation of H at the first observed survival time t1.
#'
#' @param Ht1 The solver solves for an optimal value of Ht1 by equating the estimating equation to zero.
#' @param Z The observed survival time, which is the minimum of T, C and A.
#' @param nu The censoring indicator for T or C
#' @param t A fixed time point
#' @param theta Vector of parameters
#' @param X Data matrix with covariates related to T.
#' @param W Data matrix with covariates related to C.

SolveHt1 <- function(Ht1,Z,nu,t,X,W,theta){
  k = ncol(X)
  l = ncol(W)
  beta = theta[1:k]
  eta = theta[(k+1):(l+k)]
  rho = theta[(l+k+1)]
  upx = Ht1-X%*%beta;
  upy = Ht1-W%*%eta;
  dN1 = sum((Z==t)*nu)
  c = cbind(upx,upy)
  temp = 1- pnorm(upx)-pnorm(upy)+ pbivnorm(c,rho=rho)
  return(dN1+sum((Z>=t)*log(temp)))
}


#' @title Search function
#'
#' @description Function to indicate position of t in observed survival time
#'
#' @param t fixed time t
#' @param T1 distinct observed survival time


SearchIndicate = function(t,T1){
  i = 1;
  m = length(T1);

  while(T1[i+1]<=t){
    i = i+1;
    if (i == m) break;
  }
  return(i)
}

#' @title Distance between vectors
#' @description This function computes distance between two vectors based on L2-norm
#' @param a First vector
#' @param b Second vector

Distance = function(b,a){   # L2-norm of difference
  x = b-a
  n = length(x)
  l2norm = sqrt(sum(x^2)/n)
  return(l2norm)
}

#' @title Compute bivariate survival probability
#' @description This function calculates a bivariate survival probability based on multivariate normal distribution.
#' @param lx The first lower bound of integration
#' @param ly The second lower bound
#' @param rho  Association parameter
#' @import MASS
#' @import mvtnorm

Bvprob = function(lx,ly,rho) { # compute bivariate prob.
  cor = diag(2)
  cor[1,2] = rho
  cor[2,1] = rho
  bvprob = pmvnorm(lower=c(lx,ly), upper=c(Inf,Inf), mean=c(0,0), corr = cor,algorithm = "GenzBretz");
  return(bvprob)
}

#' @title Change H to long format
#'
#' @description
#'  Change a nonparametric transformation function to long format
#'
#' @param Z Observed survival time, which is the minimum of T, C and A, where A is the administrative censoring time.
#' @param T1 Distinct observed survival time
#' @param H Nonparametric transformation function estimate

Longfun = function(Z,T1,H){
  n = length(Z)
  Hlong = rep(0,n)

  for(i in 1:n){
    j = SearchIndicate(Z[i],T1)
    Hlong[i] = H[j]
  }
return(Hlong)
}


#' @title Score equations of finite parameters
#' @description This function computes the score vectors  and the Jacobean matrix for finite model parameters.
#'
#' @inheritParams SolveH
#' @param H The estimated non-parametric transformation function for a given value of theta
#' @importFrom stats pnorm  dnorm qnorm sd

ScoreEqn = function(theta,resData,X,W,H){
  k = ncol(X)
  l = ncol(W)
  u = k+l
  beta = theta[1:k]
  eta = theta[(k+1):u]
  rho = theta[(u+1)]
  Z = resData$Z
  d1 = resData$d1
  d2 = resData$d2
  nu = d1+d2
  T1 = c(-Inf,Z[nu==1])
  T1 = unique(T1)               # distinct observed survival time
  n = length(Z);
  nu = d1+d2;

  VBeta = rep(0,k); VEta = rep(0,l);   VRho = 0;
  V = diag(1,(u+1));

  for(i in 1:n){

    up1 = H[SearchIndicate(Z[i],T1)]-X[i,]%*%beta;
    up2 = H[SearchIndicate(Z[i],T1)]-W[i,]%*%eta;

    if (up1>4) up1 = Inf
    if (up2>4) up2 = Inf
    if (up1<=-4.1) up1 = -Inf
    if (up2<=-4.1) up2 = -Inf
    up21 = (up2-up1*rho)/(1-rho^2)^0.5
    up12 = (up1-up2*rho)/(1-rho^2)^0.5
    if (up1==-Inf & up2==-Inf){
      up21 = 0
      up12 = 0
    }
    if (up1==Inf & up2==Inf){
      up21 = 0
      up12 = 0
    }
    if (up1>-Inf & up2>-Inf&up1<Inf&up2<Inf)
    {
      if(d1[i]==1 & d2[i] == 0)
      {
        SB = c(up1)*c(X[i,])-c(dnorm(up21)/(1-pnorm(up21)))*rho/(1-rho^2)^0.5*c(X[i,]);
        SE = c(1/(1-pnorm(up21))*dnorm(up21))*c(W[i,])/(1-rho^2)^0.5;
        SR = 1/(1-pnorm(up21))*(-dnorm(up21))*(rho/(1-rho^2)^1.5*up2-up1/(1-rho^2)^0.5-rho^2/(1-rho^2)^1.5*up1);
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
      if(d1[i]==0 & d2[i] == 1)
      {
        SB = c(1/(1-pnorm(up12))*dnorm(up12))*c(X[i,])/(1-rho^2)^0.5
        SE = c(up2)*W[i,]-c(dnorm(up12)/(1-pnorm(up12)))*rho/(1-rho^2)^0.5*W[i,]
        SR = 1/(1-pnorm(up12))*(-dnorm(up12))*(rho/(1-rho^2)^1.5*up1-up2/(1-rho^2)^0.5-rho^2/(1-rho^2)^1.5*up2);
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
      if(d1[i]==0 & d2[i] == 0)
      {
        temp = Bvprob(up1,up2,rho)[1];
        SB = c(dnorm(up1))*X[i,]*c(1-pnorm(up21))/temp;
        SE = c(dnorm(up2))*W[i,]*c(1-pnorm(up12))/temp;
        SR = 1/(1-rho^2)^0.5*dnorm(up1)*dnorm(up21)/temp;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
    }
    if (up1>-Inf & up2 == -Inf)
    {
      if(d1[i]==1 & d2[i] == 0)
      {
        SB = c(up1)*X[i,];
        SE = rep(0,k);
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
      if(d1[i]== 0 & d2[i] == 1)
      {
        SB = rep(0,k);
        SE = rep(0,k);
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
      if(d1[i]== 0 & d2[i] == 0)
      {
        SB = c(dnorm(up1))*X[i,]/c(1-pnorm(up1))
        SE = rep(0,k);
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
    }
    if (up1 == -Inf & up2>-Inf)
    {
      if(d1[i]== 1 & d2[i] == 0)
      {
        SB = rep(0,k);
        SE = rep(0,k);
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
      if(d1[i]== 0 & d2[i] == 1)
      {
        SB = rep(0,k);
        SE = c(up2)*W[i,];
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
      if (d1[i] == 0 & d2[i] == 0)
      {
        SB = rep(0,k);
        SE = c(dnorm(up2))*W[i,]/c(1-pnorm(up2))
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
    }
    if(up1 == -Inf & up2 == -Inf)
    {
      SB = rep(0,k);
      SE = rep(0,k);
      SR = 0;
      SR = 0;
      VBeta = VBeta+SB;
      VEta = VEta+SE;
      VRho = VRho+SR;
      S = c(SB,SE,SR);
      V = V + S%*%t(S);
    }
    if (up1 == Inf & up2>-4& up2<4)
    {
      if(d1[i] == 1 & d2[i] == 0)
      {
        SB = rep(0,k);
        SE = rep(0,k);
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
      if(d1[i]== 0 & d2[i] == 1)
      {
        SB = rep(0,k);
        SE = c(up2)*W[i,]
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
      if (d1[i] ==0 & d2[i] == 0)
      {
        SB = rep(0,k);
        SE = rep(0,k)
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
    }
    if (up1 == Inf & up2 == Inf)
    {
      if(d1[i] == 1 & d2[i] == 0)
      {
        SB = rep(0,k);
        SE = rep(0,k);
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
      if(d1[i]== 0 & d2[i] == 1)
      {
        SB = rep(0,k);
        SE = rep(0,k);
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
      if (d1[i] ==0 & d2[i] == 0)
      {
        SB = rep(0,k);
        SE = rep(0,k)
        SR = 0;
        VBeta = VBeta+SB;
        VEta = VEta+SE;
        VRho = VRho+SR;
        S = c(SB,SE,SR);
        V = V + S%*%t(S);
      }
    }
  }
  out = list(Vbeta = VBeta,Veta = VEta, Vrho = VRho,V1 = V)
  return(out)
}


#' @title  Estimate finite parameters based on score equations
#' @description This function estimates the model parameters
#'
#' @inheritParams SolveH
#' @param H The estimated non-parametric transformation function for a given value of theta.
#' @param eps Convergence error.


SolveScore = function(theta,resData,X,W,H, eps = 1e-3){   # Estimate model parameters

  PEst = ScoreEqn(theta,resData,X,W,H)

  VBeta = PEst$Vbeta
  VEta = PEst$Veta
  VRho = PEst$Vrho
  S = c(VBeta,VEta,VRho)
  V = PEst$V1
  b = theta+solve(V)%*%S
  m = length(b)
  if(b[m]>0.98) b[m] = 0.98

  step = 0
  a = theta
  while (Distance(b,a)>eps){
    a = b
    step = step+1
    PEst = ScoreEqn(a,resData,X,W,H)
    VBeta = PEst$Vbeta
    VEta = PEst$Veta
    VRho = PEst$Vrho
    S = c(VBeta,VEta,VRho)
    V = PEst$V1
    b = a+solve(V)%*%S
    if(b[m]>0.98) b[m] = 0.98    # Perfect correlation can lead to error
    if (step>=10)                # Max ten steps
      break;
  }
  return(b)
}



#' @title Fit a semiparametric transformation model for dependent censoring
#'
#' @description This function allows to estimate the dependency parameter along all other model parameters. First, estimates a non-parametric transformation function, and
#' then at the second stage it estimates other model parameters assuming that the non-parametric function is known. The details for
#' implementing the dependent censoring methodology can be found in Deresa and Van Keilegom (2021).
#'
#' @references Deresa, N. and Van Keilegom, I. (2021). On semiparametric modelling, estimation and inference for survival data subject to dependent censoring, Biometrika, 108, 965–979.
#'
#' @param start Initial values for the finite dimensional parameters. If \code{start} is NULL, the initial values will be obtained
#' by fitting an Accelerated failure time models.
#' @param resData Data matrix with three columns;  Z = the observed survival time, d1 = the censoring indicator of T
#' and  d2 =  the censoring indicator of C.
#' @param X Data matrix with covariates related to T
#' @param W Data matrix with covariates related to C
#' @param eps Convergence error. This is set by the user in such away that the desired convergence is met; the default is \code{eps = 1e-3}.
#' @param n.iter Number of iterations; the default is \code{n.iter = 20}. The larger the number of iterations, the longer the computational time.
#' @param bootstrap A boolean indicating whether to compute bootstrap standard errors for making inferences.
#' @param n.boot Number of bootstrap samples to use in the estimation of bootstrap standard errors if \code{bootstrap = TRUE}. The default is n.boot = 50. But, higher
#' values  of \code{n.boot} are recommended for obtaining good estimates of bootstrap standard errors.
#' @importFrom stats pnorm  qnorm sd
#' @importFrom survival survreg Surv
#' @importFrom MASS mvrnorm
#' @import pbivnorm nleqslv SemiPar.depCens
#'
#' @return This function returns a fit of a semiparametric transformation model; parameter estimates, estimate of the non-parametric transformation function, bootstrap standard
#' errors for finite-dimensional parameters, the nonparametric cumulative hazard function, etc.
#'
#' @examples
#' \donttest{
#' # Toy data example to illustrate implementation
#' n = 300
#' beta = c(0.5, 1); eta = c(1,1.5); rho = 0.70
#' sigma = matrix(c(1,rho,rho,1),ncol=2)
#' err = MASS::mvrnorm(n, mu = c(0,0) , Sigma=sigma)
#' err1 = err[,1]; err2 = err[,2]
#' x1 = rbinom(n,1,0.5); x2 = runif(n,-1,1)
#' X = matrix(c(x1,x2),ncol=2,nrow=n); W = X   # data matrix
#' T1 = X%*%beta+err1
#' C =  W%*%eta+err2
#' T1 = exp(T1); C = exp(C)
#' A = runif(n,0,8); Y = pmin(T1,C,A)
#' d1 = as.numeric(Y==T1)
#' d2 = as.numeric(Y==C)
#' resData = data.frame("Z" = Y,"d1" = d1, "d2" = d2)   # should be data frame
#' colnames(X) = c("X1", "X2")
#' colnames(W) = c("W1","W2")
#'
#' #  Bootstrap is false by default
#' output = NonParTrans(resData = resData, X = X, W = W, n.iter = 2)
#' output$parameterEstimates
#'
#' }
#
#' @export

NonParTrans = function(resData, X, W, start = NULL, n.iter = 15,  bootstrap = FALSE, n.boot = 50, eps = 1e-3){
  X = as.matrix(X)
  W = as.matrix(W)

  #  Verify column names of input dataframe resData
  if (!all(c("Z", "d1", "d2") %in% colnames(resData)))
    stop("Z, d1 and d2 arguments must be column names of resData")

  id <- match(c("Z", "d1", "d2"), colnames(resData))
  colnames(resData)[id] <- c("Z", "d1", "d2")

  # Make  sure that the data is ordered
  resData= resData[order(resData$Z),]
  Z = resData$Z
  d1 = resData$d1
  d2 = resData$d2
  X = X[order(Z),]
  W = W[order(Z),]

  if(is.null(start)){
    fit1 = survreg(Surv(Z,d1)~X, dist ="lognormal")
    fit2 = survreg(Surv(Z,d2)~W, dist = "lognormal")
    start = c(fit1$coefficients[-1],fit2$coefficients[-1],0.3)
  }
  l = length(start)
  if(start[l] < -1 | start[l]>1) start[l] = 0.3

  TransHat = SolveH(start, resData, X, W)
  T1 = TransHat[,1]
  H = TransHat[,2]
  parhat = SolveScore(start, resData,X,W, H, eps = eps)

  flag = 0;
  b = parhat
  a = start
  m = length(b)
  if(b[m]>0.99){
    b[m] = 0.98
  }else if (b[m] <= -0.99){
    b[m] = -0.98
    }else b = b
  while (Distance(b,a)>eps){
    a = b
    TransHat = SolveH(a, resData, X, W)
    H = TransHat[,2]
    parhat = SolveScore(a, resData, X, W, H, eps = eps)
    b = parhat
    if(b[m]>0.99){
      b[m] = 0.98
    }else if (b[m]<= -0.99){
      b[m] = -0.98
    }else b[m] = b[m]

    flag = flag+1
    if (flag>n.iter)
    {
      flag=0;
      break;
    }
  }
  TransHat = SolveH(parhat, resData, X, W)

  # Nonparametric bootstrap for making inference

  paramsBootstrap = list("init" = parhat,"resData" = resData, "X" = X, "W" = W,"n.boot" = n.boot, "n.iter" = n.iter, "eps" = eps)

  fitNonTrans <- NULL
  if(bootstrap) {                                                          # Obtain bootstrap standard error
    fitNonTrans <- do.call(boot.nonparTrans, paramsBootstrap)
    stError <- fitNonTrans$parEst.std.error
    Z <- parhat/stError
    pvalue <- 2*(1-pnorm(abs(Z)))
    out.res <- cbind(parhat,stError,pvalue)
    rownames(out.res) <- c(colnames(X), colnames(W),"rho")
    colnames(out.res) <- c("coef.", "BootSE", "Pvalue")
    output <- list("parameterEstimates" = out.res, "bootstrap" = bootstrap, "dimX" = ncol(X), "dimW" = ncol(W),"observedTime" = T1,"H" = H)
  }else{
    rownames(parhat) <- c(colnames(X), colnames(W),"rho")
    output <- c(list("parameterEstimates" = parhat, "bootstrap" = bootstrap, "dimX" = ncol(X), "dimW" = ncol(W),"observedTime" = T1,"H" = H),fitNonTrans)
  }
  class(output) <- append(class(output), "fitNonPar")                              # dependent censoring fit object

  return(output)

}


#' @title Nonparametric bootstrap approach for a Semiparametric transformation model under dependent censpring

#' @description This function estimates the bootstrap standard errors for the finite-dimensional model parameters and for the non-parametric transformation
#' function. Parallel computing using foreach has been used to speed up the estimation of standard errors.
#'
#'
#' @param init Initial values for the finite dimensional parameters obtained from the fit of \code{\link{NonParTrans}}
#' @param resData Data matrix with three columns;  Z = the observed survival time, d1 = the censoring indicator of T
#' and  d2 =  the censoring indicator of C.
#' @param X Data matrix with covariates related to T
#' @param W Data matrix with covariates related to C.
#' @param eps Convergence error. This is set by the user in such away that the desired convergence is met; the default is \code{eps = 1e-3}
#' @param n.iter Number of iterations; the default is \code{n.iter = 15}. The larger the number of iterations, the longer the computational time.
#' @param n.boot Number of bootstraps to use in the estimation of bootstrap standard errors.
#' @importFrom stats pnorm  qnorm sd
#' @importFrom survival survreg Surv
#' @import pbivnorm MASS mvtnorm nleqslv foreach parallel
#'
#' @return Bootstrap standard errors for parameter estimates and for estimated cumulative hazard function.


boot.nonparTrans = function(init,resData,X,W,n.boot, n.iter, eps){
  B = n.boot                                               # number of bootstrap samples
  n.cores <- parallel::detectCores() - 1

  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)

  boot = foreach(b = 1:B, .packages= c('survival', 'pbivnorm','nleqslv','MASS','mvtnorm'), .export = c("SolveH","SolveScore","Bvprob","Distance","ScoreEqn","Longfun" ,"SearchIndicate","SolveHt1")) %dopar% {
    samp1 = sample(length(resData$Z),replace = TRUE)
    resData_b = resData[samp1,]
    resData_b = resData_b[order(resData_b$Z),]
    Zb = resData_b$Z
    Xb = X[samp1,]
    Wb = W[samp1,]
    Xb = Xb[order(Zb),]
    Wb = Wb[order(Zb),]
    nu = resData_b$d1+resData_b$d2
    Tb1 = c(-Inf,Zb[nu==1])
    Tb1 = unique(Tb1)       # distinct observed survival time

    TransHat = SolveH(init, resData_b, Xb, Wb)
    T1 = TransHat[,1]
    H = TransHat[,2]
    parhat = SolveScore(init, resData_b,Xb,Wb, H, eps)

    flag = 0;
    b = parhat
    a = init
    m = length(b)
    if(b[m]>0.99){
      b[m] = 0.98
    }else if (b[m] <= -0.99){
      b[m] = -0.98
    }else b = b

    while (Distance(b,a)>eps){
      a = b
      TransHat = SolveH(a, resData_b, Xb, Wb)
      H = TransHat[,2]
      parhat = SolveScore(a, resData_b, Xb, Wb, H, eps)
      b = parhat
      if(b[m]>0.99){
        b[m] = 0.98
      }else if (b[m]<= -0.99){
        b[m] = -0.98
      }else b[m] = b[m]

      flag = flag+1
      if (flag>n.iter)
      {
        flag=0;
        break;
      }
    }
    longB = Longfun(Zb,Tb1,H)
    longB[longB == -Inf] = -5
    list("beta.star" = b, "H.star" = longB)
  }
  parallel::stopCluster(cl = my.cluster)

  beta.star = t(sapply(boot, function(x) x$beta.star))
  H.star =  t(sapply(boot, function(x) x$H.star))
  Bootse = apply(beta.star,2,sd)
  Hse = apply(H.star,2,sd)
  bootR = list("parEst.std.error" = Bootse, "H.std.error" = Hse)
  return(bootR)
}
