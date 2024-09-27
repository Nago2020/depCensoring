#' @title Log-likelihood function for the independence copula.
#' @description This likelihood function is maximized to estimate the model parameters under the independence copula.
#' @param para Estimated parameter values/initial values.
#' @param Y Follow-up time.
#' @param Delta Censoring indicator.
#' @param Dist.T The distribution to  be used for the survival time T. This argument can take one of the values from \code{c("lnorm", "weibull", "llogis")}.
#' @param Dist.C The distribution to  be used for the censoring time C. This argument can take one of the values from \code{c("lnorm", "weibull", "llogis")}.
#' @import rvinecopulib rafalib stats
#' @return Maximized log-likelihood value.

loglike.indep.unconstrained<-function(para,Y,Delta,Dist.T,Dist.C){
  par.T=para[1:2]
  par.C=para[3:4]
  n=length(Y)

  if(Dist.T=="lnorm"){
    FT<- plnorm(Y,meanlog=par.T[1],sdlog=exp(par.T[2]))
    fT<- (1/(sqrt(2*pi)*exp(par.T[2])*Y))*exp(-((log(Y)-par.T[1])^2/(2*(exp(par.T[2]))^2)))
  } else if(Dist.T=="weibull"){
    FT<-1-exp(-(Y/exp(par.T[2]))^(exp(par.T[1])))
    fT<-(exp(par.T[1])/exp(par.T[2]))*((Y/exp(par.T[2]))^(exp(par.T[1])-1))*(exp(-(Y/exp(par.T[2]))^(exp(par.T[1]))))
  } else if(Dist.T=="llogis"){
    FT<-1/(1+exp(-(log(Y)-par.T[1])/exp(par.T[2])))
    fT<-(1/(Y*exp(par.T[2])))*exp((log(Y)-par.T[1])/exp(par.T[2]))*(1+exp((log(Y)-par.T[1])/exp(par.T[2])))^(-2)
  }

  if(Dist.C=="lnorm"){
    FC<-plnorm(Y,meanlog=par.C[1],sdlog=exp(par.C[2]))
    fC<-(1/(sqrt(2*pi)*exp(par.C[2])*Y))*exp(-((log(Y)-par.C[1])^2/(2*(exp(par.C[2]))^2)))
  } else if(Dist.C=="weibull"){
    FC<-1-exp(-(Y/exp(par.C[2]))^(exp(par.C[1])))
    fC<-(exp(par.C[1])/exp(par.C[2]))*((Y/exp(par.C[2]))^(exp(par.C[1])-1))*(exp(-(Y/exp(par.C[2]))^(exp(par.C[1]))))
  } else if(Dist.C=="llogis"){
    FC<-1/(1+exp(-(log(Y)-par.C[1])/exp(par.C[2])))
    fC<-(1/(Y*exp(par.C[2])))*exp((log(Y)-par.C[1])/exp(par.C[2]))*(1+exp((log(Y)-par.C[1])/exp(par.C[2])))^(-2)
  }

  l1<-sum(log(fT[Delta==1])+log(1-FC[Delta==1]))
  l0<-sum(log(fC[Delta==0])+log(1-FT[Delta==0]))
  l1+l0
}

#' @title Log-likelihood function for the Frank copula.
#' @description This likelihood function is maximized to estimate the model parameters under the Frank copula.
#' @param para Estimated parameter values/initial values.
#' @param Y Follow-up time.
#' @param Delta Censoring indicator.
#' @param Dist.T The distribution to  be used for the survival time T. This argument can take one of the values from \code{c("lnorm", "weibull", "llogis")}.
#' @param Dist.C The distribution to  be used for the censoring time C. This argument can take one of the values from \code{c("lnorm", "weibull", "llogis")}.
#' @import rvinecopulib rafalib stats
#' @return Maximized log-likelihood value.

loglike.frank.unconstrained<-function(para,Y,Delta,Dist.T,Dist.C){
  par.T=para[1:2]
  par.C=para[3:4]
  par.z=para[5]
  par.tau<-tanh(par.z)
  par.cop<-sign(par.tau) * ktau_to_par("frank", max(1e-20, abs(par.tau)))
  n=length(Y)

  if(Dist.T=="lnorm"){
    FT<-plnorm(Y,meanlog=par.T[1],sdlog=exp(par.T[2]))
    fT<- (1/(sqrt(2*pi)*exp(par.T[2])*Y))*exp(-((log(Y)-par.T[1])^2/(2*(exp(par.T[2]))^2)))
  } else if(Dist.T=="weibull"){
    FT<-1-exp(-(Y/exp(par.T[2]))^(exp(par.T[1])))
    fT<-(exp(par.T[1])/exp(par.T[2]))*((Y/exp(par.T[2]))^(exp(par.T[1])-1))*(exp(-(Y/exp(par.T[2]))^(exp(par.T[1]))))
  } else if(Dist.T=="llogis"){
    FT<-1/(1+exp(-(log(Y)-par.T[1])/exp(par.T[2])))
    fT<-(1/(Y*exp(par.T[2])))*exp((log(Y)-par.T[1])/exp(par.T[2]))*(1+exp((log(Y)-par.T[1])/exp(par.T[2])))^(-2)
  }

  if(Dist.C=="lnorm"){
    FC<-plnorm(Y,meanlog=par.C[1],sdlog=exp(par.C[2]))
    fC<-(1/(sqrt(2*pi)*exp(par.C[2])*Y))*exp(-((log(Y)-par.C[1])^2/(2*(exp(par.C[2]))^2)))
  } else if(Dist.C=="weibull"){
    FC<-1-exp(-(Y/exp(par.C[2]))^(exp(par.C[1])))
    fC<-(exp(par.C[1])/exp(par.C[2]))*((Y/exp(par.C[2]))^(exp(par.C[1])-1))*(exp(-(Y/exp(par.C[2]))^(exp(par.C[1]))))
  } else if(Dist.C=="llogis"){
    FC<-1/(1+exp(-(log(Y)-par.C[1])/exp(par.C[2])))
    fC<-(1/(Y*exp(par.C[2])))*exp((log(Y)-par.C[1])/exp(par.C[2]))*(1+exp((log(Y)-par.C[1])/exp(par.C[2])))^(-2)
  }

  frank_cop<-bicop_dist("frank",0,par.cop)
  hC.condT<-hbicop(cbind(FT,FC), 1,frank_cop)
  hT.condC<-hbicop(cbind(FT,FC), 2,frank_cop)
  l1<-sum(log(fT[Delta==1])+log(1-hC.condT[Delta==1]))
  l0<-sum(log(fC[Delta==0])+log(1-hT.condC[Delta==0]))
  l1+l0
}

#' @title Log-likelihood function for the Gumbel copula.
#' @description This likelihood function is maximized to estimate the model parameters under the Gumbel copula.
#' @param para Estimated parameter values/initial values.
#' @param Y Follow-up time.
#' @param Delta Censoring indicator.
#' @param Dist.T The distribution to  be used for the survival time T. This argument can take one of the values from \code{c("lnorm", "weibull")} and has to be the same as Dist.C.
#' @param Dist.C The distribution to  be used for the censoring time C. This argument can take one of the values from \code{c("lnorm", "weibull")} and has to be the same as Dist.T.
#' @import rvinecopulib rafalib stats
#' @return Maximized log-likelihood value.

loglike.gumbel.unconstrained<-function(para,Y,Delta,Dist.T,Dist.C){
  par.T=para[1:2]
  par.C=para[3:4]
  par.z=min(c(para[5],700))
  par.tau<-exp(par.z)/(1+exp(par.z))
  par.cop<-sign(par.tau) * ktau_to_par("gumbel", max(1e-9, abs(par.tau)))
  par.cop<-min(50,par.cop)
  par.cop<-max(1,par.cop)
  n=length(Y)

  if(Dist.T=="lnorm"){
    FT<-plnorm(Y,meanlog=par.T[1],sdlog=exp(par.T[2]))
    fT<- (1/(sqrt(2*pi)*exp(par.T[2])*Y))*exp(-((log(Y)-par.T[1])^2/(2*(exp(par.T[2]))^2)))
  } else if(Dist.T=="weibull"){
    FT<-1-exp(-(Y/exp(par.T[2]))^(exp(par.T[1])))
    fT<-(exp(par.T[1])/exp(par.T[2]))*((Y/exp(par.T[2]))^(exp(par.T[1])-1))*(exp(-(Y/exp(par.T[2]))^(exp(par.T[1]))))
  }

  if(Dist.C=="lnorm"){
    FC<-plnorm(Y,meanlog=par.C[1],sdlog=exp(par.C[2]))
    fC<-(1/(sqrt(2*pi)*exp(par.C[2])*Y))*exp(-((log(Y)-par.C[1])^2/(2*(exp(par.C[2]))^2)))
  } else if(Dist.C=="weibull"){
    FC<-1-exp(-(Y/exp(par.C[2]))^(exp(par.C[1])))
    fC<-(exp(par.C[1])/exp(par.C[2]))*((Y/exp(par.C[2]))^(exp(par.C[1])-1))*(exp(-(Y/exp(par.C[2]))^(exp(par.C[1]))))
  }

  gumbel_cop<-bicop_dist("gumbel",0,par.cop)
  hC.condT<-hbicop(cbind(FT,FC), 1,gumbel_cop)
  hT.condC<-hbicop(cbind(FT,FC), 2,gumbel_cop)
  l1<-sum(log(fT[Delta==1])+log(1-hC.condT[Delta==1]))
  l0<-sum(log(fC[Delta==0])+log(1-hT.condC[Delta==0]))
  l1+l0
}

#' @title Log-likelihood function for the Clayton copula.
#' @description This likelihood function is maximized to estimate the model parameters under the Clayton copula.
#' @param para Estimated parameter values/initial values.
#' @param Y Follow-up time.
#' @param Delta Censoring indicator.
#' @param Dist.T The distribution to  be used for the survival time T. This argument can take one of the values from \code{c("lnorm", "weibull")} and has to be the same as Dist.C.
#' @param Dist.C The distribution to  be used for the censoring time C. This argument can take one of the values from \code{c("lnorm", "weibull")} and has to be the same as Dist.T.
#' @import rvinecopulib rafalib stats
#' @return Maximized log-likelihood value.

loglike.clayton.unconstrained<-function(para,Y,Delta,Dist.T,Dist.C){
  par.T=para[1:2]
  par.C=para[3:4]
  par.z=min(c(para[5],700))
  par.tau<-exp(par.z)/(1+exp(par.z))
  par.cop<-sign(par.tau) * ktau_to_par("clayton", max(1e-9, abs(par.tau)))
  par.cop<-min(28,par.cop)
  par.cop<-max(1e-10,par.cop)
  n=length(Y)

  if(Dist.T=="lnorm"){
    FT<-plnorm(Y,meanlog=par.T[1],sdlog=exp(par.T[2]))
    fT<- (1/(sqrt(2*pi)*exp(par.T[2])*Y))*exp(-((log(Y)-par.T[1])^2/(2*(exp(par.T[2]))^2)))
  } else if(Dist.T=="weibull"){
    FT<-1-exp(-(Y/exp(par.T[2]))^(exp(par.T[1])))
    fT<-(exp(par.T[1])/exp(par.T[2]))*((Y/exp(par.T[2]))^(exp(par.T[1])-1))*(exp(-(Y/exp(par.T[2]))^(exp(par.T[1]))))
  }

  if(Dist.C=="lnorm"){
    FC<-plnorm(Y,meanlog=par.C[1],sdlog=exp(par.C[2]))
    fC<-(1/(sqrt(2*pi)*exp(par.C[2])*Y))*exp(-((log(Y)-par.C[1])^2/(2*(exp(par.C[2]))^2)))
  } else if(Dist.C=="weibull"){
    FC<-1-exp(-(Y/exp(par.C[2]))^(exp(par.C[1])))
    fC<-(exp(par.C[1])/exp(par.C[2]))*((Y/exp(par.C[2]))^(exp(par.C[1])-1))*(exp(-(Y/exp(par.C[2]))^(exp(par.C[1]))))
  }

  clayton_cop<-bicop_dist("clayton",0,par.cop)
  hC.condT<-hbicop(cbind(FT,FC), 1,clayton_cop)
  hT.condC<-hbicop(cbind(FT,FC), 2,clayton_cop)
  l1<-sum(log(fT[Delta==1])+log(1-hC.condT[Delta==1]))
  l0<-sum(log(fC[Delta==0])+log(1-hT.condC[Delta==0]))
  l1+l0
}

#' @title Log-likelihood function for the Gaussian copula.
#' @description This likelihood function is maximized to estimate the model parameters under the Gaussian copula.
#' @param para Estimated parameter values/initial values.
#' @param Y Follow-up time.
#' @param Delta Censoring indicator.
#' @param Dist.T The distribution to  be used for the survival time T. This argument can only the value \code{"lnorm"}.
#' @param Dist.C The distribution to  be used for the censoring time C. This argument can only the value \code{"lnorm"}.
#' @import rvinecopulib rafalib stats
#' @return Maximized log-likelihood value.

loglike.gaussian.unconstrained<-function(para,Y,Delta,Dist.T,Dist.C){
  par.T=para[1:2]
  par.C=para[3:4]
  par.z=para[5]
  par.tau<-tanh(par.z)
  par.cop<-sign(par.tau) * ktau_to_par("gaussian", max(1e-20, abs(par.tau)))
  n=length(Y)

  if(Dist.T=="lnorm"){
    FT<-plnorm(Y,meanlog=par.T[1],sdlog=exp(par.T[2]))
    fT<- (1/(sqrt(2*pi)*exp(par.T[2])*Y))*exp(-((log(Y)-par.T[1])^2/(2*(exp(par.T[2]))^2)))
  } else if(Dist.T=="weibull"){
    FT<-1-exp(-(Y/exp(par.T[2]))^(exp(par.T[1])))
    fT<-(exp(par.T[1])/exp(par.T[2]))*((Y/exp(par.T[2]))^(exp(par.T[1])-1))*(exp(-(Y/exp(par.T[2]))^(exp(par.T[1]))))
  }

  if(Dist.C=="lnorm"){
    FC<-plnorm(Y,meanlog=par.C[1],sdlog=exp(par.C[2]))
    fC<-(1/(sqrt(2*pi)*exp(par.C[2])*Y))*exp(-((log(Y)-par.C[1])^2/(2*(exp(par.C[2]))^2)))
  } else if(Dist.C=="weibull"){
    FC<-1-exp(-(Y/exp(par.C[2]))^(exp(par.C[1])))
    fC<-(exp(par.C[1])/exp(par.C[2]))*((Y/exp(par.C[2]))^(exp(par.C[1])-1))*(exp(-(Y/exp(par.C[2]))^(exp(par.C[1]))))
  }

  gaussian_cop<-bicop_dist("gaussian",0,par.cop)
  hC.condT<-hbicop(cbind(FT,FC), 1,gaussian_cop)
  hT.condC<-hbicop(cbind(FT,FC), 2,gaussian_cop)
  l1<-sum(log(fT[Delta==1])+log(1-hC.condT[Delta==1]))
  l0<-sum(log(fC[Delta==0])+log(1-hT.condC[Delta==0]))
  l1+l0
}

#' @title Fit the dependent censoring models.
#' @description Estimates the model parameters by maximizing the log-likelihood.
#' @param Y Follow-up time.
#' @param Delta Censoring indicator.
#' @param Copula The copula family. This argument can take values from \code{c("frank","gumbel","clayton","gaussian","indep")}.
#' @param Dist.T The distribution to  be used for the survival time T. This argument can take one of the values from \code{c("lnorm", "weibull", "llogis")}.
#' @param Dist.C The distribution to  be used for the censoring time C. This argument can take one of the values from \code{c("lnorm", "weibull", "llogis")}.
#' @import rvinecopulib rafalib stats
#' @return A list containing the minimized negative log-likelihood using the independence copula model, the estimated parameter values for the model with the independence copula, the minimized negative log-likelihood using the specified copula model and the estimated parameter values for the model with the specified copula.

optimlikelihood<-function(Y,Delta,Copula,Dist.T,Dist.C){

  para.start=optim(c(1,1,1,1),loglike.indep.unconstrained,Y=Y,Delta=Delta,Dist.T=Dist.T,Dist.C=Dist.C,control=list(fnscale=-1),method="BFGS")

  if(Copula=="frank"){
    result=optim(c(para.start$par,0),loglike.frank.unconstrained,Y=Y,Delta=Delta,Dist.T=Dist.T,Dist.C=Dist.C,control=list(fnscale=-1),method="BFGS")
    return(list(para.start$value,para.start$par,result$value,result$par))
  } else if(Copula=="clayton"){
    result=optim(c(para.start$par,0),loglike.clayton.unconstrained,Y=Y,Delta=Delta,Dist.T=Dist.T,Dist.C=Dist.C,control=list(fnscale=-1),method="BFGS")
    return(list(para.start$value,para.start$par,result$value,result$par))
  } else if(Copula=="gumbel"){
    result=optim(c(para.start$par,0),loglike.gumbel.unconstrained,Y=Y,Delta=Delta,Dist.T=Dist.T,Dist.C=Dist.C,control=list(fnscale=-1),method="BFGS")
    return(list(para.start$value,para.start$par,result$value,result$par))
  } else if(Copula=="gaussian"){
    result=optim(c(para.start$par,0),loglike.gaussian.unconstrained,Y=Y,Delta=Delta,Dist.T=Dist.T,Dist.C=Dist.C,control=list(fnscale=-1),method="BFGS")
    return(list(para.start$value,para.start$par,result$value,result$par))
  } else if(Copula == "indep"){
    return(list(para.start$value,para.start$par))
  }
}

#' @title Estimation of the specified parametric dependent censoring model without covariates. Note that it is not assumed that the association parameter
#' of the copula function is known, unlike most other papers in the literature. The details for implementing the methodology can be found in Czado and Van Keilegom (2023).
#' @references Czado and Van Keilegom (2023). Dependent censoring based on parametric copulas. Biometrika, 110(3), 721-738.
#' @description Returns the estimated the model parameters.
#' @inheritParams optimlikelihood
#' @import rvinecopulib rafalib stats
#' @return A table containing the minimized negative log-likelihood using the independence copula model, the estimated parameter values for the model with the independence copula, the minimized negative log-likelihood using the specified copula model and the estimated parameter values for the model with the specified copula.
#'
#' @examples
#' \donttest{
#'tau = 0.75
#'Copula = "frank"
#'Dist.T = "weibull"
#'Dist.C = "lnorm"
#'par.T = c(2,1)
#'par.C = c(1,2)
#'n=1000
#'
#'simdata<-TCsim(tau,Copula,Dist.T,Dist.C,par.T,par.C,n)
#'
#'Y = simdata[[1]]
#'Delta = simdata[[2]]
#'
#'ParamCop(Y,Delta,Copula,Dist.T,Dist.C)
#' }
#'
#' @export

ParamCop<-function(Y,Delta,Copula,Dist.T,Dist.C){

  if(!(Copula %in% c("frank","gumbel","clayton","gaussian","indep"))){
    print("ERROR : Copula needs to be frank, gumbel, clayton, gaussian or indep")
    return(NA)
  }else if(!(Dist.T %in% c("lnorm","weibull","llogis"))){
    print("ERROR : Dist.T needs to be lnorm, weibull or llogis")
    return(NA)
  }else if(!(Dist.C %in% c("lnorm","weibull","llogis"))){
    print("ERROR : Dist.C needs to be lnorm, weibull or llogis")
    return(NA)
  }else if( (Copula == "gumbel" | Copula == "clayton" | Copula == "gaussian") & (Dist.T == "llogis" | Dist.C == "llogis")){
    print("ERROR : log-logistic marginals cannot be specified with the gumbel, clayton or gaussian copula")
    return(NA)
  }else if( (Copula == "gumbel" | Copula == "clayton" | Copula == "gaussian") & (Dist.T != Dist.C)){
    print("ERROR : marginals need to be both lnorm or both weibull")
    return(NA)
  } else {

    result = optimlikelihood(Y,Delta,Copula,Dist.T,Dist.C)

    if (Copula == "indep"){
      par.mat<-rbind(c(result[[1]],result[[2]]))
      rownames(par.mat)<-c("indep")
      colnames(par.mat)<-c("log-like","E(log(T))","sd(log(T))","E(log(C))","sd(log(C))")
      par.mat[,3]<-exp(par.mat[,3])
      par.mat[,5]<-exp(par.mat[,5])
      if (Dist.T=="weibull"){
        par.mat[,2]<-exp(par.mat[,2])
        colnames(par.mat)[2]<-c("shape(T)")
        colnames(par.mat)[3]<-c("scale(T)")
      }
      if (Dist.C=="weibull"){
        par.mat[,4]<-exp(par.mat[,4])
        colnames(par.mat)[4]<-c("shape(C)")
        colnames(par.mat)[5]<-c("scale(C)")
      }
      if (Dist.T=="llogis"){
        colnames(par.mat)[2]<-c("location(log(T))")
        colnames(par.mat)[3]<-c("scale(log(T))")
      }
      if (Dist.C=="llogis"){
        colnames(par.mat)[4]<-c("location(log(C))")
        colnames(par.mat)[5]<-c("scale(log(C))")
      }

    } else {
      par.mat<-rbind(c(result[[1]],result[[2]],0),c(result[[3]],result[[4]]))

      rownames(par.mat)<-c("indep",Copula)
      colnames(par.mat)<-c("log-like","E(log(T))","sd(log(T))","E(log(C))","sd(log(C))","tau")

      par.mat[,3]<-exp(par.mat[,3])
      par.mat[,5]<-exp(par.mat[,5])

      if (Dist.T=="weibull"){
        par.mat[,2]<-exp(par.mat[,2])
        colnames(par.mat)[2]<-c("shape(T)")
        colnames(par.mat)[3]<-c("scale(T)")
      }
      if (Dist.C=="weibull"){
        par.mat[,4]<-exp(par.mat[,4])
        colnames(par.mat)[4]<-c("shape(C)")
        colnames(par.mat)[5]<-c("scale(C)")
      }

      if (Dist.T=="llogis"){
        colnames(par.mat)[2]<-c("location(log(T))")
        colnames(par.mat)[3]<-c("scale(log(T))")
      }
      if (Dist.C=="llogis"){
        colnames(par.mat)[4]<-c("location(log(C))")
        colnames(par.mat)[5]<-c("scale(log(C))")
      }

      if (Copula=="frank" | Copula =="gaussian"){
        par.mat[2,6]<-tanh(par.mat[2,6])
      }
      if (Copula=="clayton" | Copula =="gumbel"){
        par.mat[2,6]<-exp(par.mat[2,6])/(1+exp(par.mat[2,6]))
      }
    }
    return(round(par.mat,digits=3))
  }
}


#' @title Function to simulate (Y,Delta) from the copula based model for (T,C).
#' @description Generates the follow-up time and censoring indicator according to the specified model.
#' @param tau Value of Kendall's tau for (T,C). The default value is 0.
#' @param Copula The copula family. This argument can take values from \code{c("frank","gumbel","clayton","gaussian","indep")}. The default copula model is "frank".
#' @param Dist.T Distribution of the survival time T. This argument can take one of the values from \code{c("lnorm", "weibull", "llogis")}. The default distribution is "lnorm".
#' @param Dist.C Distribution of the censoring time C. This argument can take one of the values from \code{c("lnorm", "weibull", "llogis")}. The default distribution is "lnorm".
#' @param par.T Parameter values for the distribution of T.
#' @param par.C Parameter values for the distribution of C.
#' @param n Sample size.
#' @return A list containing the generated follow-up times and censoring indicators.
#' @import rvinecopulib rafalib stats
#'
#' @examples
#' \donttest{
#' tau = 0.5
#' Copula = "gaussian"
#' Dist.T = "lnorm"
#' Dist.C = "lnorm"
#' par.T = c(1,1)
#' par.C = c(2,2)
#' n=1000
#'
#' simdata <- TCsim(tau,Copula,Dist.T,Dist.C,par.T,par.C,n)
#' Y = simdata[[1]]
#' Delta = simdata[[2]]
#'
#' hist(Y)
#' mean(Delta)
#' }
#'
#' @export

TCsim<-function(tau=0,Copula="frank",Dist.T="lnorm",Dist.C="lnorm",
                par.T=c(0,1),par.C=c(0,1),n=10000)
{

  if(Copula=="indep"){
    u.TC = cbind(runif(n),runif(n))
  } else {
    parTC<-ktau_to_par(Copula,tau)
    u.TC<-rbicop(n,Copula,0, parTC)
  }

  f.TC<-rep(0,n)
  F.TC<-rep(0,n)
  hC.condT<-rep(0,n)
  hT.condC<-rep(0,n)

  if(Dist.T=="lnorm"){
    xT<-qlnorm(u.TC[,1],meanlog=par.T[1],sdlog=par.T[2])
  }
  else if(Dist.T=="weibull"){
    xT<-qweibull(u.TC[,1],par.T[1],par.T[2])
  }
  else if(Dist.T=="llogis"){
    xT<-exp(qlogis(u.TC[,1],par.T[1],par.T[2]))
  }

  if(Dist.C=="lnorm"){
    xC<-qlnorm(u.TC[,2],meanlog=par.C[1],sdlog=par.C[2])
  }
  else if(Dist.C=="weibull"){
    xC<-qweibull(u.TC[,2],par.C[1],par.C[2])
  }
  else if(Dist.C=="llogis"){
    xC<-exp(qlogis(u.TC[,2],par.C[1],par.C[2]))
  }
  Delta<-as.numeric(xT<=xC)
  Y<-pmin(xT,xC)
  return(list(Y,Delta))
}




