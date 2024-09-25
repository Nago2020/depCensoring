
#' @title Competing risk likelihood function.
#'
#' @description This function implements the second step likelihood function of
#' the competing risk model defined in Rutten, Willems et al. (202X).
#'
#' @param n The sample size.
#' @param s The number of competing risks.
#' @param Y The observed times.
#' @param admin Boolean value indicating whether or not administrative censoring
#' should be taken into account.
#' @param cens.inds matrix of censoring indicators (each row corresponds to a
#' single observation).
#' @param M Design matrix, consisting of [intercept, exo.cov, Z, cf]. Note that
#' \code{cf} represents the multiple ways of 'handling' the endogenous covariate
#' Z, see also the documentation of 'estimate.cmprsk.R'. When there is no
#' confounding, \code{M} will be [intercept, exo.cov].
#' @param Sigma The covariance matrix.
#' @param beta.mat Matrix containing all of the covariate effects.
#' @param sigma.vct Vector of standard deviations. Should be equal to
#' \code{sqrt(diag(Sigma))}.
#' @param rho.mat The correlation matrix.
#' @param theta.vct Vector containing the parameters of the Yeo-Johnsontrans-
#' formations.
#'
#' @import mvtnorm pbivnorm OpenMx stats
#'
#' @return Evaluation of the log-likelihood function
#'

cr.lik <- function(n, s, Y, admin, cens.inds, M, Sigma, beta.mat, sigma.vct,
                   rho.mat, theta.vct) {

  # In the current implementation, Sigma should always be a positive definite
  # matrix. The if-clause is therefore superfluous.
  if (is.positive.definite(Sigma, tol = 1e-30)) {

    # Compute the (derivate of the) Yeo-Johnson transformations of Y
    transY.T <- matrix(nrow = n, ncol = s)
    DtransY.T <- matrix(nrow = n, ncol = s)
    for (i in 1:s) {
      transY.T[, i] <- YJtrans(Y, theta.vct[i])
      DtransY.T[, i] <- DYJtrans(Y, theta.vct[i])
    }

    # Precompute some values used in obtaining b_T and m_{k, j}.
    Mbeta <- M %*% beta.mat

    # Compute b_Tj/sigma_j, for each j.
    b_T <- (transY.T - Mbeta)/matrix(rep(sigma.vct, n), nrow = n, byrow = TRUE)

    # For each latent time T^j, compute the arguments to be supplied to the
    # multivariate normal distribution function.
    args.f_T <- array(Inf, dim = c(n, s, s + as.numeric(admin)))
    partial.correlations <- array(Inf, dim = c(s, s, s))
    for (j in 1:s) {

      # Compute the value of m_{k.j} for each observation
      mj <- matrix(Inf, nrow = n, ncol = s)
      for (k in 1:s) {
        if (k != j) {
          mj[, k] <- Mbeta[, k] + rho.mat[j, k] * sigma.vct[k] * b_T[, j]
        }
      }
      # Compute s_{k.j} for each observation
      sj <- rep(Inf, s)
      for (k in 1:s) {
        if (k != j) {
          sj[k] <- sigma.vct[k]*(1 - rho.mat[j, k]^2)^(1/2)
        }
      }

      # Compute rho_{kq.j}
      rhoj <- matrix(Inf, nrow = s, ncol = s)
      for (k in 1:s) {
        for (q in 1:s) {
          if ((k != j) & (q != j)) {
            rhoj[k, q] <- (rho.mat[k, q] - rho.mat[j, k]*rho.mat[j, q])/sqrt((1 - rho.mat[j, k]^2)*(1 - rho.mat[j, q]^2))
          }
        }
      }

      # Arguments to be supplied to the multivariate normal distribution
      args.f_Tj <- -(transY.T - mj)/matrix(rep(sj, n), nrow = n, byrow = TRUE)
      args.f_T[, , j] <- args.f_Tj

      # Partial correlations
      partial.correlations[, , j] <- rhoj
    }

    # Also compute the likelihood contributions when the times would correspond to
    # administrative censoring, if necessary
    if (admin) {
      args.f_Ta <- -b_T/matrix(rep(sigma.vct, n), nrow = n, byrow = TRUE)
      args.f_T[, , s + 1] <- args.f_Ta
    }

    # Matrix of arguments for which the multivariate normal distribution should
    # be evaluated.
    args.mvn <- matrix(nrow = n, ncol = s)
    for (i in 1:s) {
      args.mvn[, i] <- apply(args.f_T[, i ,]^cens.inds, 1, prod)
    }

    # Define function to be used to evaluate the multivariate normal distr.
    # function. The choice can greatly impact the computation time.
    if (s == 3) {
      mypmvnorm.cr <- function(arg.vct, cov.mat) {
        pbivnorm(arg.vct[1], arg.vct[2], rho = cov.mat[1, 2]/sqrt(cov.mat[1, 1]*cov.mat[2, 2]))
      }
    } else {
      mypmvnorm.cr <- function(arg.vct, cov.mat) {
        omxMnor(cov.mat, rep(0, length(arg.vct)), rep(-Inf, length(arg.vct)), arg.vct)
      }
    }
    mypmvnorm.ad <- function(arg.vct, cov.mat) {
      omxMnor(cov.mat, rep(0, length(arg.vct)), rep(-Inf, length(arg.vct)), arg.vct)
    }

    # Multivariate normal evaluations
    mvn.evals <- rep(Inf, n)
    cov.mat.array <- array(dim = c(s, s, s + as.numeric(admin)))
    cov.mat.array[, , 1:s] <- partial.correlations
    if (admin) {
      cov.mat.array[, , s + 1] <- rho.mat
    }

    for (row.idx in 1:nrow(args.mvn)) {
      event.type.idx <- which(cens.inds[row.idx, ] == 1)
      if (event.type.idx <= s) { # If it is a competing risks
        cov.mat <- cov.mat.array[-event.type.idx, -event.type.idx, event.type.idx]
        mvn.evals[row.idx] <- mypmvnorm.cr(args.mvn[row.idx, !is.na(args.mvn[row.idx, ])], cov.mat)
      } else {                  # If it is administrative censoring
        cov.mat <- cov.mat.array[, , event.type.idx]
        mvn.evals[row.idx] <- mypmvnorm.ad(args.mvn[row.idx, ], cov.mat)
      }
    }

    # Compute the likelihood contributions
    lik.contr <- rep(Inf, n)
    for (i in 1:n) {
      event.type.idx <- which(cens.inds[i, ] == 1)
      if (event.type.idx <= s) { # If it is a competing risks
        lik.contr[i] <- sigma.vct[event.type.idx]^(-1) * mvn.evals[i] * dnorm(b_T[i, event.type.idx]) * DtransY.T[i, event.type.idx]
      } else {                  # If it is administrative censoring
        lik.contr[i] <- mvn.evals[i]
      }
    }

    # Return the likelihood
    lik.contr <- pmax(lik.contr, 1e-100)
    Logn <- sum(log(lik.contr))
    return(-Logn)
  } else {
    return(2^{60}-1)                # if not positive definite --> return a very large value
  }
}



#' @title First step log-likelihood function for Z continuous
#'
#' @description This function defines the maximum likelihood used to estimate
#' the control function in the case of a continuous endogenous variable.
#'
#' @param gamma Vector of coefficients in the linear model used to estimate the
#' control function.
#' @param Z Endogenous covariate.
#' @param M Design matrix, containing an intercept, the exogenous covariates and
#' the instrumental variable.
#'
#' @return Returns the log-likelihood function corresponding to the data,
#' evaluated at the point \code{gamma}.
#'

LikGamma1 <- function(gamma, Z, M) {

  # Correctly format the variables
  W <- as.matrix(M)
  gamma <- as.matrix(gamma)

  # Vector of squared errors in the linear regression model
  tot <- (Z - W %*% gamma)^2

  # Prevent numerical issues
  p1 <- pmax(tot, 1e-100)

  # Sum of (log of) squared errors to be minimized.
  Logn <- sum(log(p1));

  # Return the results
  return(Logn)
}



#' @title First step log-likelihood function for Z binary.
#'
#' @description This function defines the maximum likelihood used to estimate
#' the control function in the case of a binary endogenous variable.
#'
#' @param gamma Vector of coefficients in the logistic model used to estimate
#' the control function.
#' @param Z Endogenous covariate.
#' @param M Design matrix, containing an intercept, the exogenous covariates and
#' the instrumental variable.
#'
#' @import stats
#'
#' @return Returns the log-likelihood function corresponding to the data,
#' evaluated at the point \code{gamma}.
#'

LikGamma2 <- function(gamma, Z, M) {

  # Correcly format the variables
  W <- as.matrix(M)
  gamma <- as.matrix(gamma)

  # Factors in likelihood function of logistic regression
  tot <- (plogis(W %*% gamma)^Z)*((1 - plogis(W %*% gamma))^(1-Z))

  # Prevent numerical issues
  p1 <- pmax(tot, 1e-100)

  # Log-likelihood function to be maximized when finding gamma
  Logn <- sum(log(p1))

  # In order to be consistent with 'LikGamma1.R' (which should be minimized to
  # find gamma), we return the negative of the log-likelihood, which should also
  # be minimized to find gamma.
  return(-Logn)
}



#' @title Second step log-likelihood function.
#'
#' @description This function defines the log-likelihood used to estimate
#' the second step in the competing risks extension of the model described in
#' Rutten, Willems et al. (20XX).
#'
#' @param par Vector of all second step model parameters, consisting of the
#' regression parameters, variance-covariance matrix elements and transformation
#' parameters.
#' @param data Data frame resulting from the 'uniformize.data.R' function.
#' @param admin Boolean value indicating whether the data contains
#' administrative censoring.
#' @param conf Boolean value indicating whether the data contains confounding
#' and hence indicating the presence of Z and W.
#' @param cf "Control function" to be used. This can either be the (i) estimated
#' control function, (ii) the true control function, (iii) the instrumental
#' variable, or (iv) nothing (\code{cf = NULL}). Option (ii) is used when
#' comparing the two-step estimator to the oracle estimator, and option (iii) is
#' used to compare the two-step estimator with the naive estimator.
#'
#' @import mvtnorm pbivnorm OpenMx
#'
#' @return Log-likelihood evaluation of the second step.
#'

LikF.cmprsk <- function(par, data, admin, conf, cf) {

  #### Unpack data ####

  # Observed times
  Y <- data[, "y"]

  # Censoring indicators
  cens.inds <- data[, grepl("delta[[:digit:]]+", colnames(data))]
  if (admin) {
    cens.inds <- cbind(cens.inds, data[, which(colnames(data) == "da")])
  }

  # Intercept
  intercept <- data[, "x0"]

  # Exogenous covariates
  exo.cov <- data[, grepl("x[1-9][[:digit:]]*", colnames(data)),
                  drop = FALSE]

  # Endogenous covariate
  if (conf) {
    Z <- data[, "z"]
  } else {
    Z <- NULL
    cf <- NULL
  }

  # Design matrix: [intercept, exogenous covariates, Z, cf]
  M <- as.matrix(cbind(intercept, as.matrix(exo.cov), Z, cf))

  #### Unpack parameters ####

  # Define some useful variables
  k <- ncol(M)                                             # Nbr. parameters
  s <- ifelse(admin, ncol(cens.inds) - 1, ncol(cens.inds)) # Nbr. comp. risks
  n <- nrow(M)                                             # Nbr. observations

  # Extract parameters for beta
  beta.vct <- par[1:(k*s)]
  beta.mat <- matrix(nrow = k, ncol = s)
  for (i in 1:s) {
    beta.mat[,i] <- beta.vct[((i - 1)*k + 1):(i*k)]
  }

  # Extract parameters for the variance-covariance matrix.
  # Structure of sd-vector:
  # [1:s]       -> standerd deviations
  # [(s+1):end] -> correlation; stored as
  #                [rho_{12}, rho_{13}, ..., rho_{1d}, rho_{23}, rho_{24}, ...]
  sd <- par[(k*s + 1):(length(par) - s)]

  # Construct correlation and covariance matrix
  rho.mat <- matrix(1, nrow = s, ncol = s)
  rho.mat[lower.tri(rho.mat, diag = FALSE)] <- sd[(s + 1):length(sd)]
  rho.mat[upper.tri(rho.mat, diag = FALSE)] <- t(rho.mat)[upper.tri(rho.mat, diag = FALSE)]
  Sigma <- rho.mat * outer(sd[1:s], sd[1:s], "*")
  sigma.vct <- sd[1:s]

  # Extract parameters for transformation function
  theta.vct <- par[(length(par) - s + 1):(length(par))]

  #### Compute log-likelihood evaluation ####

  cr.lik(n, s, Y, admin, cens.inds, M, Sigma, beta.mat, sigma.vct, rho.mat, theta.vct)
}



#' @title Wrapper implementing likelihood function using Cholesky factorization.
#'
#' @description This function parametrizes the covariance matrix using its
#' Cholesky decomposition, so that optimization of the likelihood can be done
#' based on this parametrization, and positive-definiteness of the covariance
#' matrix is guaranteed at each step of the optimization algorithm.
#'
#' @param par.chol Vector of all second step model parameters, consisting of the
#' regression parameters, Cholesky decomposition of the variance-covariance
#' matrix elements and transformation parameters.
#' @param data Data frame resulting from the 'uniformize.data.R' function.
#' @param admin Boolean value indicating whether the data contains
#' administrative censoring.
#' @param conf Boolean value indicating whether the data contains confounding
#' and hence indicating the presence of Z and W.
#' @param cf "Control function" to be used. This can either be the (i) estimated
#' control function, (ii) the true control function, (iii) the instrumental
#' variable, or (iv) nothing (\code{cf = NULL}). Option (ii) is used when
#' comparing the two-step estimator to the oracle estimator, and option (iii) is
#' used to compare the two-step estimator with the naive estimator.
#' @param eps Minimum value for the diagonal elements in the covariance matrix.
#' Default is \code{eps = 0.001}.
#'
#' @import mvtnorm pbivnorm OpenMx
#'
#' @return Log-likelihood evaluation of the second step.
#'

likF.cmprsk.Cholesky <- function(par.chol, data, admin, conf, cf, eps = 0.001) {

  #### Unpack data ####

  # Censoring indicators
  cens.inds <- data[, grepl("delta[[:digit:]]+", colnames(data))]
  if (admin) {
    cens.inds <- cbind(cens.inds, data[, which(colnames(data) == "da")])
  }

  # Intercept
  intercept <- data[, "x0"]

  # Exogenous covariates
  exo.cov <- data[, grepl("x[1-9][[:digit:]]*", colnames(data)),
                  drop = FALSE]

  # Endogenous covariate and control function
  if (conf) {
    Z <- data[, "z"]
  } else {
    Z <- NULL
    cf <- NULL
  }

  # Design matrix: [intercept, exogenous covariates, Z, estimated cf]
  M <- as.matrix(cbind(intercept, as.matrix(exo.cov), Z, cf))

  #### Transform Cholesky factorization parameters ####

  # Define some useful variables
  k <- ncol(M)                                             # Nbr. parameters
  s <- ifelse(admin, ncol(cens.inds) - 1, ncol(cens.inds)) # Nbr. comp. risks

  # Extract parameters for the variance-covariance matrix. Note that in this
  # case, the parameters of the Cholesky factorization are provided.
  sd.chol <- par.chol[(k*s + 1):(length(par.chol) - s)]

  # Reconstruct upper diagonal matrix of Cholesky factorization
  chol.U <- matrix(0, nrow = s, ncol = s)
  chol.U[lower.tri(chol.U, diag = TRUE)] <- sd.chol

  # Construct variance-covariance matrix
  Sigma <- chol.U %*% t(chol.U)
  diag(Sigma) <- diag(Sigma) + eps

  # Correlation matrix
  rho.mat <- Sigma / outer(sqrt(diag(Sigma)), sqrt(diag(Sigma)), FUN = "*")

  # Get the variances and covariances
  sd <- rep(0, length(sd.chol))
  sd[1:s] <- sqrt(diag(Sigma))
  sd[(s+1):length(sd)] <- t(rho.mat)[lower.tri(rho.mat)]

  # Contruct vector to be supplied to likF.cmprsk.R
  par <- par.chol
  par[(k*s + 1):(length(par) - s)] <- sd

  #### Run regular likelihood function ####

  LikF.cmprsk(par, data, admin, conf, cf)
}



#' @title Second likelihood function needed to fit the independence model in the
#' second step of the estimation procedure.
#'
#' @description This function defines the log-likelihood used in estimating
#' the second step in the competing risks extension of the model described in
#' Rutten, Willems et al. (20XX). The results of this function will serve as
#' starting values for subsequent optimizations (LikI.comprsk.R and
#' LikF.cmprsk.R)
#'
#' @param par Vector of all second step model parameters, consisting of the
#' regression parameters, variance-covariance matrix elements and transformation
#' parameters.
#' @param data Data frame resulting from the 'uniformize.data.R' function.
#' @param admin Boolean value indicating whether the data contains
#' administrative censoring.
#' @param conf Boolean value indicating whether the data contains confounding
#' and hence indicating the presence of Z and W
#' @param cf "Control function" to be used. This can either be the (i) estimated
#' control function, (ii) the true control function, (iii) the instrumental
#' variable, or (iv) nothing (\code{cf = NULL}). Option (ii) is used when
#' comparing the two-step estimator to the oracle estimator, and option (iii) is
#' used to compare the two-step estimator with the naive estimator.
#'
#' @import stats
#'
#' @return Starting values for subsequent optimization function used in the
#' second step of the estimation procedure.
#'

LikI.bis <- function(par, data, admin, conf, cf) {

  #### Unpack data ####

  # Observed times
  Y <- data[, "y"]

  # Censoring indicators
  cens.inds <- data[, grepl("delta[[:digit:]]+", colnames(data))]
  cens.inds <- cbind(cens.inds, data[, which(colnames(data) == "da")])

  # Intercept
  intercept <- data[, "x0"]

  # Exogenous covariates
  exo.cov <- data[, grepl("x[1-9][[:digit:]]*", colnames(data)),
                  drop = FALSE]

  # Endogenous covariate
  if (conf) {
    Z <- data[, "z"]
  } else {
    Z <- NULL
    cf <- NULL
  }

  # Design matrix: [intercept, exogenous covariates, Z, cf]
  M <- as.matrix(cbind(intercept, as.matrix(exo.cov), Z, cf))

  #### Unpack parameters ####

  # Define some useful variables
  k <- ncol(M)
  s <- ifelse(admin, ncol(cens.inds) - 1, ncol(cens.inds))
  n <- nrow(M)

  # Extract parameters for beta
  beta.vct <- par[1:(k*s)]
  beta.mat <- matrix(nrow = k, ncol = s)
  for (i in 1:s) {
    beta.mat[,i] <- beta.vct[((i - 1)*k + 1):(i*k)]
  }

  # Extract parameters for the variance-covariance matrix
  sigma.vct <- par[(k*s + 1):(k*s + s)]

  # Construct variance-covariance matrix (in this function, all correlations are
  # assumed to equal zero)
  Sigma <- diag(sigma.vct[1:s]^2, nrow = s)

  # Extract parameters for transformation function
  theta.vct <- par[(k*s + s + 1):(length(par))]

  #### Compute the likelihood function evaluation ####

  # For computational efficiency, we de not use the function 'cr.lik.contr.R' to
  # obtain the likelihood contributions. Computational speed-up can be achieved
  # by exploiting the independence assumption, which allows to evaluate the
  # multivariate normal distribution in a more efficient way.

  # Compute the (derivate of the) Yeo-Johnson transformations of Y
  transY.T <- matrix(nrow = n, ncol = s)
  DtransY.T <- matrix(nrow = n, ncol = s)
  for (i in 1:s) {
    transY.T[, i] <- YJtrans(Y, theta.vct[i])
    DtransY.T[, i] <- DYJtrans(Y, theta.vct[i])
  }

  # Compute b_Tj, for each j.
  b_T <- matrix(nrow = n, ncol = s)
  for (i in 1:s) {
    b_T[, i] <- (transY.T[,i] - (M %*% beta.mat[,i]))/sigma.vct[i]
  }

  # Precompute the values m_{k, j}. Note that these values are equal for all j
  # because of the independence (second equation on page 22).
  m <- M %*% beta.mat

  # For each latent time T^j, compute the likelihood contribution
  f_T <- matrix(nrow = n, ncol = s)
  for (j in 1:s) {
    mj <- m
    mj[, j] <- Inf
    args <- (transY.T - mj)/matrix(rep(sigma.vct, n), nrow = n, byrow = TRUE)
    f_T[,j] <- 1/sigma.vct[j] * dnorm(b_T[, j]) * DtransY.T[, j] *
      apply(pnorm(-args), 1, prod)
  }

  # Also compute the likelihood contributions when the times would correspond to
  # administrative censoring
  if (admin) {
    args <- b_T/matrix(rep(sigma.vct, n), nrow = n, byrow = TRUE)
    f_Ta <- apply(pnorm(-args), 1, prod)
    f_T <- cbind(f_T, f_Ta)
  }

  tot <- apply(f_T^cens.inds, 1, prod)
  p1 <- pmax(tot, 1e-100)
  Logn <- sum(log(p1))
  return(-Logn)
}



#' @title Second step log-likelihood function under independence assumption.
#'
#' @description This function defines the log-likelihood used to estimate
#' the second step in the competing risks extension assuming independence of
#' some of the competing risks in the model described in Rutten, Willems et al.
#' (20XX).
#'
#' @param par Vector of all second step model parameters, consisting of the
#' regression parameters, variance-covariance matrix elements and transformation
#' parameters.
#' @param data Data frame resulting from the 'uniformize.data.R' function.
#' @param eoi.indicator.names Vector of names of the censoring indicator columns
#' pertaining to events of interest. Events of interest will be modeled allowing
#' dependence between them, whereas all censoring events (corresponding to
#' indicator columns not listed in \code{eoi.indicator.names}) will be treated
#' as independent of every other event.
#' @param admin Boolean value indicating whether the data contains
#' administrative censoring.
#' @param conf Boolean value indicating whether the data contains confounding
#' and hence indicating the presence of Z and W
#' @param cf "Control function" to be used. This can either be the (i) estimated
#' control function, (ii) the true control function, (iii) the instrumental
#' variable, or (iv) nothing (\code{cf = NULL}). Option (ii) is used when
#' comparing the two-step estimator to the oracle estimator, and option (iii) is
#' used to compare the two-step estimator with the naive estimator.
#'
#' @import mvtnorm pbivnorm OpenMx
#'
#' @return Log-likelihood evaluation for the second step in the esimation
#' procedure.
#'

LikI.cmprsk <- function(par, data, eoi.indicator.names, admin, conf, cf) {

  #### Unpack data ####

  # Observed times
  Y <- data[, "y"]

  # Censoring indicators
  cens.inds <- data[, grepl("delta[[:digit:]]+", colnames(data))]
  if (admin) {
    cens.inds <- cbind(cens.inds, data[, which(colnames(data) == "da")])
  }

  # Intercept
  intercept <- data[, "x0"]

  # Exogenous covariates
  exo.cov <- data[, grepl("x[1-9][[:digit:]]*", colnames(data)),
                  drop = FALSE]

  # Endogenous covariate and instrumental variable
  if (conf) {
    Z <- data[, "z"]
  } else {
    Z <- NULL
    cf <- NULL
  }

  # Design matrix: [intercept, exogenous covariates, Z, estimated cf]
  M <- as.matrix(cbind(intercept, as.matrix(exo.cov), Z, cf))

  #### Unpack parameters ####

  # Define some useful variables
  k <- ncol(M)                                             # Nbr. parameters
  s <- ifelse(admin, ncol(cens.inds) - 1, ncol(cens.inds)) # Nbr. comp. risks
  n <- nrow(M)                                             # Nbr. observations

  # Extract parameters for beta
  beta.vct <- par[1:(k*s)]
  beta.mat <- matrix(nrow = k, ncol = s)
  for (i in 1:s) {
    beta.mat[,i] <- beta.vct[((i - 1)*k + 1):(i*k)]
  }

  # Extract parameters for the variance-covariance matrix. Note that in this
  # case, the correlations between the competing risks and the censoring is zero.
  # Structure of sd-vector:
  # [1:s]       -> standerd deviations
  # [(s+1):end] -> correlation; stored as
  #                [rho_{ij}, rho_{ik}, ..., rho_{iq} rho_{jk}, ...],
  #                where eoi.inds = c(i, j, k, ..., q)
  sd <- par[(k*s + 1):(length(par) - s)]

  # Quick compatibility check
  eoi.inds <- as.numeric(gsub("delta", "", eoi.indicator.names))
  if (length(sd) != s + choose(length(eoi.inds), 2)) {
    stop("Supplied parameter vector is not consistent with necessary model parameters.")
  }

  # Construct correlation and covariance matrix
  rho.mat <- matrix(0, nrow = s, ncol = s)
  diag(rho.mat) <- 1
  sd.idx <- s + 1
  for (i in eoi.inds) {
    for (j in eoi.inds) {
      if (i < j) {
        rho.mat[i, j] <- sd[sd.idx]
        rho.mat[j, i] <- sd[sd.idx]
        sd.idx <- sd.idx + 1
      }
    }
  }
  Sigma <- rho.mat * outer(sd[1:s], sd[1:s], FUN = "*")
  sigma.vct <- sd[1:s]

  # Extract parameters for transformation function
  theta.vct <- par[(length(par) - s + 1):(length(par))]

  #### Compute log-likelihood evaluation ####

  cr.lik(n, s, Y, admin, cens.inds, M, Sigma, beta.mat, sigma.vct, rho.mat, theta.vct)

}


#' @title Wrapper implementing likelihood function assuming independence between
#' competing risks and censoring using Cholesky factorization.
#'
#' @description This function does the same as LikI.cmprsk (in fact, it even
#' calls said function), but it parametrizes the covariance matrix using its
#' Cholesky decomposition in order to guarantee positive definiteness. This
#' function is never used, might not work and could be deleted.
#'
#' @param par.chol Vector of all second step model parameters, consisting of the
#' regression parameters, Cholesky decomposition of the variance-covariance
#' matrix elements and transformation parameters.
#' @param data Data frame resulting from the 'uniformize.data.R' function.
#' @param eoi.indicator.names Vector of names of the censoring indicator columns
#' pertaining to events of interest. Events of interest will be modeled allowing
#' dependence between them, whereas all censoring events (corresponding to
#' indicator columns not listed in \code{eoi.indicator.names}) will be treated
#' as independent of every other event.
#' @param admin Boolean value indicating whether the data contains
#' administrative censoring.
#' @param conf Boolean value indicating whether the data contains confounding
#' and hence indicating the presence of Z and W.
#' @param cf "Control function" to be used. This can either be the (i) estimated
#' control function, (ii) the true control function, (iii) the instrumental
#' variable, or (iv) nothing (\code{cf = NULL}). Option (ii) is used when
#' comparing the two-step estimator to the oracle estimator, and option (iii) is
#' used to compare the two-step estimator with the naive estimator.
#' @param eps Minimum value for the diagonal elements in the covariance matrix.
#' Default is \code{eps = 0.001}.
#'
#' @import mvtnorm pbivnorm OpenMx
#'
#' @return Log-likelihood evaluation for the second step in the estimation
#' procedure.


LikI.cmprsk.Cholesky <- function(par.chol, data, eoi.indicator.names, admin,
                                 conf, cf, eps = 0.001) {

  ## Unpack data ##

  # Censoring indicators
  cens.inds <- data[, grepl("delta[[:digit:]]+", colnames(data))]
  if (admin) {
    cens.inds <- cbind(cens.inds, data[, which(colnames(data) == "da")])
  }

  # Intercept
  intercept <- data[, "x0"]

  # Exogenous covariates
  exo.cov <- data[, grepl("x[1-9][[:digit:]]*", colnames(data)),
                  drop = FALSE]

  # Endogenous covariate and control function
  if (conf) {
    Z <- data[, "z"]
  } else {
    Z <- NULL
    cf <- NULL
  }

  # Design matrix: [intercept, exogenous covariates, Z, estimated cf]
  M <- as.matrix(cbind(intercept, as.matrix(exo.cov), Z, cf))

  #### Transform Cholesky factorization parameters ####

  # Define some useful variables
  k <- ncol(M)                                             # Nbr. parameters
  s <- ifelse(admin, ncol(cens.inds) - 1, ncol(cens.inds)) # Nbr. comp. risks

  # Extract parameters for the variance-covariance matrix. Note the particular
  # form of the covariance matrix, which is (assuming for simplicity of
  # representation that the first s1 events are the events of interest):
  #
  #      s1        s2
  #  ---------  ---------
  # [s x ... x  0 0 ... 0]
  #           ...
  # [x s ... x  0 0 ... 0]
  #
  # [0 0 ... 0  s 0 ... 0]
  # [0 0 ... 0  0 s ... 0]
  #           ...
  # [0 0 ... 0  0 0 ... s]
  #
  # where 'x' represents a possibly non-zero element and 's' represents a
  # positive element

  # Let s1 equal the number of events of interest. The vector sd.chol has the
  # following structure:
  #
  # [1:(choolse(s1, 2) + s1]: Parameters of the Cholesky decomposition of the
  #                           sub matrix of the covariance matrix pertaining to
  #                           the events of interest.
  # [choose(s1, 2) + s1:end]: Standard deviations pertaining to the events that
  #                           are to be modeled independently from every other
  #                           event.
  sd.chol <- par.chol[(k*s + 1):(length(par.chol) - s)]

  # Extract parameters for the events of interest
  s1 <- length(eoi.indicator.names)
  sd1.chol <- sd.chol[1:(choose(s1, 2) + s1)]
  chol1.U <- matrix(0, nrow = s1, ncol = s1)
  chol1.U[lower.tri(chol1.U, diag = TRUE)] <- sd1.chol
  Sigma1 <- chol1.U %*% t(chol1.U) + diag(eps, nrow = s1)
  sd1 <- sqrt(diag(Sigma1))
  rho.mat1 <- Sigma1 / outer(sd1, sd1)
  cor1 <- t(rho.mat1)[lower.tri(rho.mat1)]

  # Extract parameters for the independent events
  sd2 <- NULL
  if (s1 < s) {
    sd2 <- sd.chol[(choose(s1, 2) + s1 + 1):length(sd.chol)]
  }

  sd <- c(sd1, sd2, cor1)

  # Construct vector to be supplied to likI.cmprsk.R
  par <- par.chol
  par[(k*s + 1):(length(par.chol) - s)] <- sd

  #### Run regular likelihood function ####

  LikI.cmprsk(par, data, eoi.indicator.names, admin, conf, cf)
}



#' @title Full likelihood (including estimation of control function).
#'
#' @description This function defines the 'full' likelihood of the model.
#' Specifically, it includes the estimation of the control function in the
#' computation of the likelihood. This function is used in the estimation of the
#' variance of the estimates (variance.cmprsk.R).
#'
#' @param parhatG The full parameter vector.
#' @param data Data frame.
#' @param eoi.indicator.names Vector of names of the censoring indicator columns
#' pertaining to events of interest. Events of interest will be modeled allowing
#' dependence between them, whereas all censoring events (corresponding to
#' indicator columns not listed in \code{eoi.indicator.names}) will be treated
#' as independent of every other event. If \code{eoi.indicator.names == NULL},
#' all events will be modelled dependently.
#' @param admin Boolean value indicating whether the data contains
#' administrative censoring.
#' @param conf Boolean value indicating whether the data contains confounding
#' and hence indicating the presence of Z and W.
#' @param Zbin Boolean value indicating whether the confounding variable is
#' binary.
#' @param inst Type of instrumental function to be used.
#'
#' @import mvtnorm pbivnorm OpenMx
#'
#' @return Full model log-likelihood evaluation.
#'

likIFG.cmprsk.Cholesky <- function(parhatG, data, eoi.indicator.names, admin,
                                   conf, Zbin, inst) {

  #### Precondition checks ####

  if (!conf) {
    stop("Only implemented when there is confounding in the data set.")
  }
  if (inst != "cf") {
    stop("Only implemented when the control function should be estimated.")
  }

  #### Extract parameters ####

  # Censoring indicators
  cens.inds <- data[, grepl("delta[[:digit:]]+", colnames(data))]
  if (admin) {
    cens.inds <- cbind(cens.inds, data[, which(colnames(data) == "da")])
  }

  # Design matrix of exogenous covariates
  intercept <- data[, "x0"]
  exo.cov <- data[, grepl("x[1-9][[:digit:]]*", colnames(data)), drop = FALSE]
  W <- data[, "w"]
  XandW <- as.matrix(cbind(intercept, exo.cov, W))

  # Endogenous covariate
  Z <- data[, "z"]

  # Set useful variables and extract parameter vector
  k <- ncol(XandW) + 1                                     # Nbr. parameters
  s <- ifelse(admin, ncol(cens.inds) - 1, ncol(cens.inds)) # Nbr. comp. risks
  par.chol <- parhatG[1:(length(parhatG) - k + 1)]
  gammaest <- parhatG[(length(parhatG) - k + 2):(length(parhatG))]

  #### Estimate control function ####

  cf <- estimate.cf(XandW, Z, Zbin, gammaest)$cf

  #### Return full likelihood ####

  if (length(eoi.indicator.names) == s) {
    likF.cmprsk.Cholesky(par.chol, data, admin, conf, cf)
  } else {
    LikI.cmprsk.Cholesky(par.chol, data, eoi.indicator.names, admin, conf, cf)
  }

}



#' @title Estimate the control function
#'
#' @description This function estimates the control function for the endogenous
#' variable based on the provided covariates. This function is called inside
#' estimate.cmprsk.R.
#'
#' @param XandW Design matrix of exogenous covariates.
#' @param Z Endogenous covariate.
#' @param Zbin Boolean value indicating whether endogenous covariate is binary.
#' @param gammaest Vector of pre-estimated parameter vector. If \code{NULL},
#' this function will first estimate \code{gammaest}. Default value is
#' \code{gammaest = NULL}.
#' @import nloptr
#' @importFrom stats lm
#'
#' @return List containing the vector of values for the control function and
#' the regression parameters of the first step.
#'

estimate.cf <- function(XandW, Z, Zbin, gammaest = NULL) {

  # Define useful parameter
  parlgamma <- ncol(XandW)

  # Fit appropriate model for Z and obtain control function
  if (Zbin) {
    if (is.null(gammaest)) {
      gammaest <- nloptr(x0 = rep(0, parlgamma),
                         eval_f = LikGamma2,
                         Z = Z,
                         M = XandW,
                         lb = c(rep(-Inf, parlgamma)),
                         ub = c(rep(Inf, parlgamma)),
                         eval_g_ineq = NULL,
                         opts = list(algorithm = "NLOPT_LN_BOBYQA",
                                     "ftol_abs" = 1.0e-30,
                                     "maxeval" = 100000,
                                     "xtol_abs" = rep(1.0e-30))
      )$solution
    }
    cf <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
  } else {
    if (is.null(gammaest)) {
      gammaest <- lm(Z~XandW - 1)$coefficients
    }
    cf <- Z-(XandW%*%gammaest)
  }

  # Return the results
  list(cf = cf, gammaest = gammaest)
}



#' @title Estimate the competing risks model of Rutten, Willems et al. (20XX).
#'
#' @description This function estimates the parameters in the competing risks
#' model described in Rutten, Willems et al. (20XX). Note that this model
#' extends the model of Crommen, Beyhum, Van Keilegom (2024) and as such, this
#' function also implements their methodology.
#'
#' @param data A data frame, adhering to the following formatting rules:
#' \itemize{
#'    \item The first column, named \code{"y"}, contains the observed times.
#'    \item The next columns, named \code{"delta1"}, \code{delta2}, etc. contain
#'    the indicators for each of the competing risks.
#'    \item The next column, named \code{da}, contains the censoring indicator
#'    (independent censoring).
#'    \item The next column should be a column of all ones (representing the
#'    intercept), names \code{x0}.
#'    \item The subsequent columns should contain the values of the covariates,
#'    named \code{x1}, \code{x2}, etc.
#'    \item When applicable, the next column should contain the values of the
#'    endogenous variable. This column should be named \code{z}.
#'    \item When \code{z} is provided and an instrument for \code{z} is
#'    available, the next column, named \code{w}, should contain the values for
#'    the instrument.
#' }
#' @param admin Boolean value indicating whether the data contains
#' administrative censoring.
#' @param conf Boolean value indicating whether the data contains confounding
#' and hence indicating the presence of \code{z} and, possibly, \code{w}.
#' @param eoi.indicator.names Vector of names of the censoring indicator columns
#' pertaining to events of interest. Events of interest will be modeled allowing
#' dependence between them, whereas all censoring events (corresponding to
#' indicator columns not listed in \code{eoi.indicator.names}) will be treated
#' as independent of every other event. If \code{eoi.indicator.names == NULL},
#' all events will be modeled dependently.
#' @param Zbin Indicator value indicating whether (\code{Zbin = TRUE}) or not
#' \code{Zbin = FALSE} the endogenous covariate is binary. Default is
#' \code{Zbin = NULL}, corresponding to the case when \code{conf == FALSE}.
#' @param inst Variable encoding which approach should be used for dealing with
#' the confounding. \code{inst = "cf"} indicates that the control function
#' approach should be used. \code{inst = "W"} indicates that the instrumental
#' variable should be used 'as is'. \code{inst = "None"} indicates that Z will
#' be treated as an exogenous covariate. Finally, when \code{inst = "oracle"},
#' this function will access the argument \code{realV} and use it as the values
#' for the control function. Default is \code{inst = "cf"}.
#' @param realV Vector of numerics with length equal to the number of rows in
#' \code{data}. Used to provide the true values of the instrumental function
#' to the estimation procedure.
#' @param eps Value that will be added to the diagonal of the covariance matrix
#' during estimation in order to ensure strictly positive variances.
#'
#' @import OpenMx matrixcalc nloptr stringr numDeriv pbivnorm mvtnorm
#'
#' @return A list of parameter estimates in the second stage of the estimation
#' algorithm (hence omitting the estimates for the control function), as well
#' as an estimate of their variance and confidence intervals.
#'
#' @examples
#' \donttest{
#'
#' n <- 200
#'
#' # Set parameters
#' gamma <- c(1, 2, 1.5, -1)
#' theta <- c(0.5, 1.5)
#' eta1 <- c(1, -1, 2, -1.5, 0.5)
#' eta2 <- c(0.5, 1, 1, 3, 0)
#'
#' # Generate exogenous covariates
#' x0 <- rep(1, n)
#' x1 <- rnorm(n)
#' x2 <- rbinom(n, 1, 0.5)
#'
#' # Generate confounder and instrument
#' w <- rnorm(n)
#' V <- rnorm(n, 0, 2)
#' z <- cbind(x0, x1, x2, w) %*% gamma + V
#' realV <- z - (cbind(x0, x1, x2, w) %*% gamma)
#'
#' # Generate event times
#' err <- mvrnorm(n, mu = c(0, 0), Sigma =
#' matrix(c(3, 1, 1, 2), nrow = 2, byrow = TRUE))
#' bn <- cbind(x0, x1, x2, z, realV) %*% cbind(eta1, eta2) + err
#' Lambda_T1 <- bn[,1]; Lambda_T2 <- bn[,2]
#' T1 <- IYJtrans(Lambda_T1, theta[1]); T2 <- IYJtrans(Lambda_T2, theta[2])
#'
#' # Generate adminstrative censoring time
#' C <- runif(n, 0, 40)
#'
#' # Create observed data set
#' y <- pmin(T1, T2, C)
#' delta1 <- as.numeric(T1 == y)
#' delta2 <- as.numeric(T2 == y)
#' da <- as.numeric(C == y)
#' data <- data.frame(cbind(y, delta1, delta2, da, x0, x1, x2, z, w))
#' colnames(data) <- c("y", "delta1", "delta2", "da", "x0", "x1", "x2", "z", "w")
#'
#' # Estimate the model
#' admin <- TRUE                # There is administrative censoring in the data.
#' conf <- TRUE                 # There is confounding in the data (z)
#' eoi.indicator.names <- NULL  # We will not impose that T1 and T2 are independent
#' Zbin <- FALSE                # The confounding variable z is not binary
#' inst <- "cf"                 # Use the control function approach
#' # Since we don't use the oracle estimator, this argument is ignored anyway
#' realV <- NULL
#' estimate.cmprsk(data, admin, conf, eoi.indicator.names, Zbin, inst, realV)
#'
#' }
#'
#' @export
#'

estimate.cmprsk <- function(data, admin, conf,
                            eoi.indicator.names = NULL, Zbin = NULL,
                            inst = "cf", realV = NULL, eps = 0.001) {

  #### Conformity checks ####

  if (conf) {
    if (!is.logical(Zbin)) {
      stop("When conf = TRUE, a valid value for Zbin should be supplied.")
    }
    if (!("z" %in% colnames(data))) {
      stop("If there is confounding, the confounding variable should be names 'z'.")
    }
    if (!((inst %in% c("cf", "W", "None")) | is.vector(inst))) {
      stop("When 'conf' is TRUE, inst should be 'cf', 'W', 'None' or a vector (see documentation).")
    }
  } else {
    if (!is.null(Zbin)) {
      stop("Argument mismatch: Zbin is supplied (not NULL), yet conf = FALSE.")
    }
  }
  if ((inst == "oracle") & (is.null(realV))) {
    stop("When inst == oracle, a vector of values to be used as the control function should be provided through realV.")
  }

  # Indicator whether to use Cholesky decomposition for estimating the covariance
  # matrix. (Function will currently throw an error if use.chol == FALSE).
  use.chol <- TRUE

  #### Unpack data ####

  # Observed times
  Y <- data[, "y"]

  # Censoring indicators
  cens.inds <- data[, grepl("delta[[:digit:]]+", colnames(data))]
  cens.inds <- cbind(cens.inds, data[, which(colnames(data) == "da")])

  # Intercept
  intercept <- data[, "x0"]

  # Exogenous covariates
  exo.cov <- data[, grepl("x[1-9][[:digit:]]*", colnames(data)),
                  drop = FALSE]

  # Endogenous covariate and Instrumental variable
  if (conf) {
    Z <- data[, "z"]
    W <- data[, "w"]
  }

  # Indices of the events of interest
  if (is.null(eoi.indicator.names)) {
    eoi.indicator.names <- colnames(data)[grepl("delta[[:digit:]]", colnames(data))]
  }
  eoi.inds <- as.numeric(gsub("delta", "", eoi.indicator.names))

  # Design matrix without endogenous covariate
  XandW <- as.matrix(cbind(intercept, exo.cov))
  if (conf) {
    XandW <- cbind(XandW, W)
  }

  #### Step 1: estimate control function if necessary ####

  if (conf & (inst == "cf")) {
    est.cf.out <- estimate.cf(XandW, Z, Zbin)
    cf <- est.cf.out$cf
    gammaest <- est.cf.out$gammaest
  } else if (conf & (inst == "w")) {
    cf <- W
  } else if (conf & (inst == "none")) {
    cf <- NULL
  } else if (conf & (inst == "oracle")) {
    cf <- realV
  } else {
    cf <- NULL
  }

  #### Step 2: estimate model parameters ####

  # Define useful parameters
  n.trans <- sum(grepl("delta[1-9][[:digit:]]*", colnames(data)))
  if (conf & !is.null(cf)) {
    # When there is confounding, we have [intercept, exo.cov, Z, cf]
    totparl <- (1 + ncol(exo.cov) + 2) * n.trans
  } else if (conf & is.null(cf)) {
    # When there is confounding but no cf, we have [intercept, exo.cov, Z]
    totparl <- (1 + ncol(exo.cov) + 1) * n.trans
  } else {
    # When there is no confounding, we have [intercept exo.cov]
    totparl <- (1 + ncol(exo.cov)) * n.trans
  }

  sd.lb <- 1e-05
  sd.ub <- Inf
  cor.lb <- -0.99
  cor.ub <- 0.99
  chol.lb <- -Inf
  chol.ub <- Inf
  transfo.param.lb <- 0
  transfo.param.ub <- 2

  # Initial estimate for the parameters assuming independence between competing
  # risks.
  init.covariate.effects <- rep(0, totparl)
  init.sds <- rep(1, n.trans)
  init.transfo.pararms <- rep(1, n.trans)
  init.par <- c(init.covariate.effects, init.sds, init.transfo.pararms)

  lb.vec <- c(rep(-Inf, totparl), rep(sd.lb, n.trans), rep(transfo.param.lb, n.trans))
  ub.vec <- c(rep(Inf, totparl), rep(sd.ub, n.trans), rep(transfo.param.ub, n.trans))

  parhat.init <- nloptr(x0 = init.par,
                        eval_f = LikI.bis,
                        data = data,
                        admin = admin,
                        conf = conf,
                        cf = cf,
                        lb = lb.vec,
                        ub = ub.vec,
                        eval_g_ineq = NULL,
                        opts = list(algorithm = "NLOPT_LN_BOBYQA",
                                    "ftol_abs" = 1.0e-30,
                                    "maxeval" = 100000,
                                    "xtol_abs" = rep(1.0e-30))
                        )$solution

  ## Estimate the model parameters, possibly under the assumption of solely
  ## independent censoring if indicated. Use the pre-estimates of the model that
  ## assumes independent competing risks.

  # Set pre-estimated values as initial values
  pre.covariate.effects <- parhat.init[1:totparl]
  pre.sds <- parhat.init[(totparl + 1):(totparl + n.trans)]
  pre.transfo.params <- parhat.init[(totparl + n.trans + 1):length(parhat.init)]

  # Estimate the model
  if (!use.chol) {

    # Prevent this option from being used, as it is deprecated
    stop("Cholesky factorization should be used in estimating the model.")

    # Add initial value for the correlations between the competing risks
    pre.correlations <- rep(0, choose(length(eoi.indicator.names), 2))
    init <- c(pre.covariate.effects, pre.sds, pre.correlations, pre.transfo.params)

    lb.vec <- c(rep(-Inf, totparl), rep(sd.lb, n.trans), rep(cor.lb, choose(n.trans, 2)),
                rep(transfo.param.lb, n.trans))
    ub.vec <- c(rep(Inf, totparl), rep(sd.ub, n.trans), rep(cor.ub, choose(n.trans, 2)),
                rep(transfo.param.ub, n.trans))

    # Obtain the model parameters
    parhat1 <- nloptr(x0 = init,
                      eval_f = LikI.cmprsk,
                      data = data,
                      eoi.indicator.names = eoi.indicator.names,
                      admin = admin,
                      conf = conf,
                      cf = cf,
                      lb = lb.vec,
                      ub = ub.vec,
                      eval_g_ineq = NULL,
                      opts = list(algorithm = "NLOPT_LN_BOBYQA",
                                  "ftol_abs" = 1.0e-30,
                                  "maxeval" = 100000,
                                  "xtol_abs" = 1.0e-30)
                      )$solution

    } else {

    ## Add initial value for the correlations between the competing risks
    pre.cov.mat <- diag(pre.sds[eoi.inds]^2, nrow = length(eoi.inds),
                        ncol = length(eoi.inds))
    pre.chol1 <- chol(pre.cov.mat)[lower.tri(pre.cov.mat, diag = TRUE)]
    pre.sd2 <- pre.sds[setdiff(1:n.trans, eoi.inds)]
    init <- c(pre.covariate.effects, pre.chol1, pre.sd2, pre.transfo.params)

    s1 <- length(eoi.inds)
    s2 <- n.trans - s1
    lb.vec <- c(rep(-Inf, totparl), rep(chol.lb, choose(s1, 2) + s1),
                rep(sd.lb, s2), rep(transfo.param.lb, n.trans))
    ub.vec <- c(rep(Inf, totparl), rep(chol.ub, choose(s1, 2) + s1),
                rep(sd.ub, s2), rep(transfo.param.ub, n.trans))

    ## Obtain the model parameters
    parhatc <- nloptr(x0 = init,
                      eval_f = LikI.cmprsk.Cholesky,
                      data = data,
                      eoi.indicator.names = eoi.indicator.names,
                      admin = admin,
                      conf = conf,
                      cf = cf,
                      eps = eps,
                      lb = lb.vec,
                      ub = ub.vec,
                      eval_g_ineq = NULL,
                      opts = list(algorithm = "NLOPT_LN_BOBYQA",
                                  "ftol_abs" = 1.0e-30,
                                  "maxeval" = 100000,
                                  "xtol_abs" = 1.0e-30)
    )$solution

    ## Obtain variance of estimates

    var.out <- variance.cmprsk(parhatc, gammaest, data, admin, conf, inst, cf,
                               eoi.indicator.names, Zbin, use.chol, n.trans,
                               totparl)

    ## Backtransform the cholesky parameters

    # Extract parameters
    par.chol <- parhatc[(totparl + 1):(length(parhatc) - n.trans)]
    par.chol1 <- par.chol[1:(choose(s1, 2) + s1)]
    par.sd2 <- par.chol[(choose(s1, 2) + s1 + 1):length(par.chol)]

    # Reconstruct covariance matrix for events of interest
    mat.chol1 <- matrix(0, nrow = s1, ncol = s1)
    mat.chol1[lower.tri(mat.chol1, diag = TRUE)] <- par.chol1
    var1 <- mat.chol1 %*% t(mat.chol1)

    # Reconstruct full covariance matrix and correlation matrix
    var <- matrix(0, nrow = n.trans, ncol = n.trans)
    var[eoi.inds, eoi.inds] <- var1
    diag(var)[setdiff(1:n.trans, eoi.inds)] <- par.sd2^2
    rho.mat <- var / sqrt(outer(diag(var), diag(var), FUN = "*"))

    # Construct full parameter vector.
    parhat1 <- c(parhatc[1:totparl], diag(var),
                 rho.mat[lower.tri(rho.mat, diag = FALSE)],
                 parhatc[(length(parhatc) - n.trans + 1):length(parhatc)])

    # parhat1 is not returned as it is already computed during variance
    # estimation.

    # Store the results
    res <- var.out
    }

  #### Return the results ####

  res
}



#' @title Transform Cholesky decomposition to covariance matrix parameter element.
#'
#' @description This function transforms the parameters of the Cholesky de-
#' composition to a covariance matrix element. This function is used in
#' chol2par.R.
#'
#' @param a The row index of the covariance matrix element to be computed.
#' @param b The column index of the covariance matrix element to be computed.
#' @param par.chol1 The vector of Cholesky parameters.
#'
#' @return Specified element of the covariance matrix resulting from the
#' provided Cholesky decomposition.
#'
chol2par.elem <- function(a, b, par.chol1) {

  # Create matrix with Cholesky parameters
  p <- (-1 + sqrt(1 + 8*length(par.chol1)))/2
  chol.mat <- matrix(0, nrow = p, ncol = p)
  chol.mat[lower.tri(chol.mat, diag = TRUE)] <- par.chol1
  chol.mat[upper.tri(chol.mat, diag = TRUE)] <- t(chol.mat)[upper.tri(chol.mat, diag = TRUE)]

  # Return covariance element
  cov.elem <- 0
  for (i in 1:min(a, b)) {
    cov.elem <- cov.elem + chol.mat[a, i]*chol.mat[b, i]
  }
  cov.elem
}



#' @title Transform Cholesky decomposition to covariance matrix
#'
#' @description This function transforms the parameters of the Cholesky de-
#' composition to the covariance matrix, represented as a the row-wise con-
#' catenation of the upper-triangular elements.
#'
#' @param par.chol1 The vector of Cholesky parameters.
#'
#' @return Covariance matrix corresponding to the provided Cholesky decomposition.
chol2par <- function(par.chol1) {

  # Create matrix with Cholesky parameters
  p <- (-1 + sqrt(1 + 8*length(par.chol1)))/2
  chol.mat <- matrix(0, nrow = p, ncol = p)
  chol.mat[lower.tri(chol.mat, diag = TRUE)] <- par.chol1
  chol.mat[upper.tri(chol.mat, diag = TRUE)] <- t(chol.mat)[upper.tri(chol.mat, diag = TRUE)]

  # Define matrix of double indices corresponding to the single index rep-
  # resentation of the elements in par.chol.
  double.idx.mat <- matrix(1, nrow = 1, ncol = 2)
  for (i in 2:length(par.chol1)) {
    second <- ifelse(double.idx.mat[nrow(double.idx.mat), 2] < p,
                     double.idx.mat[nrow(double.idx.mat), 2] + 1,
                     double.idx.mat[nrow(double.idx.mat), 1] + 1)
    first <- ifelse(double.idx.mat[nrow(double.idx.mat), 2] < p,
                    double.idx.mat[nrow(double.idx.mat), 1],
                    double.idx.mat[nrow(double.idx.mat), 1] + 1)
    double.idx.mat <- rbind(double.idx.mat, c(first, second))
  }

  # Switch columns in double.idx.mat
  double.idx.mat <- cbind(double.idx.mat[, 2], double.idx.mat[, 1])

  # Compute transformation
  par.vec <- rep(NA, nrow(double.idx.mat))
  for (i in 1:nrow(double.idx.mat)) {
    par.vec[i] <- chol2par.elem(double.idx.mat[i, 1], double.idx.mat[i, 2],
                                par.chol1)
  }

  # Give appropriate row and column names
  names(par.vec) <- apply(double.idx.mat,
                          1,
                          function(row) {sprintf("(%s)", paste(row, collapse= ", "))})

  # Return the result
  par.vec
}



#' @title Derivative of transform Cholesky decomposition to covariance matrix
#' element.
#'
#' @description This function defines the derivative of the transformation
#' function that maps Cholesky parameters to a covariance matrix element. This
#' function is used in dchol2par.R.
#'
#' @param k The row index of the parameter with respect to which to take the
#' derivative.
#' @param q the column index of the parameter with respect to which to take the
#' derivative.
#' @param a The row index of the covariance matrix element to be computed.
#' @param b The column index of the covariance matrix element to be computed.
#' @param par.chol1 The vector of Cholesky parameters.
#'
#' @return Derivative of the function that transforms the cholesky parameters
#' to the specified element of the covariance matrix, evaluated at the specified
#' arguments.

dchol2par.elem <- function(k, q, a, b, par.chol1) {

  # Create matrix with Cholesky parameters
  p <- (-1 + sqrt(1 + 8*length(par.chol1)))/2
  chol.mat <- matrix(0, nrow = p, ncol = p)
  chol.mat[lower.tri(chol.mat, diag = TRUE)] <- par.chol1
  chol.mat[upper.tri(chol.mat, diag = TRUE)] <- t(chol.mat)[upper.tri(chol.mat, diag = TRUE)]

  # Return derivative of the transformation function, evaluated as par.chol1
  if (k != a & q != b) {
    der <- 0
  } else if (k == a & q == b & a == b) {
    der <- 2*chol.mat[k, q]
  } else if (k == a & k == q & a > b) {
    der <- 0
  } else if (k == a & q == b & k > q) {
    der <- chol.mat[q, q]
  } else if (k == q & a > b & b == q) {
    der <- chol.mat[a, b]
  } else if (k > q & q == b & a == b) {
    der <- 0
  } else if (k > q & q < b & a == b & k == a) {
    der <- 2*chol.mat[k, q]
  } else if (k > q & q < b & a == b & k != a) {
    der <- 0
  } else if (k > q & q == b & b < a & k != a) {
    der <- 0
  } else if (k > q & k == a & a > b & q < b) {
    der <- chol.mat[b, q]
  } else if (k > q & k == a & a > b & q > b) {
    der <- 0
  } else if (k > q & k == a & a > b & q == b) {
    der <- 0
  }
  der
}



#' @title Derivative of transform Cholesky decomposition to covariance matrix.
#'
#' @description This function defines the derivative of the transformation
#' function that maps Cholesky parameters to the full covariance matrix.
#'
#' @param par.chol1 The vector of Cholesky parameters.
#'
#' @return Derivative of the function that transforms the cholesky parameters
#' to the full covariance matrix.

dchol2par <- function(par.chol1) {

  # Create matrix with Cholesky parameters
  p <- (-1 + sqrt(1 + 8*length(par.chol1)))/2
  chol.mat <- matrix(0, nrow = p, ncol = p)
  chol.mat[lower.tri(chol.mat, diag = TRUE)] <- par.chol1
  chol.mat[upper.tri(chol.mat, diag = TRUE)] <- t(chol.mat)[upper.tri(chol.mat, diag = TRUE)]

  # Define matrix of double indices corresponding to the single index rep-
  # resentation of the elements in par.chol.
  double.idx.mat <- matrix(1, nrow = 1, ncol = 2)
  for (i in 2:length(par.chol1)) {
    second <- ifelse(double.idx.mat[nrow(double.idx.mat), 2] < p,
                     double.idx.mat[nrow(double.idx.mat), 2] + 1,
                     double.idx.mat[nrow(double.idx.mat), 1] + 1)
    first <- ifelse(double.idx.mat[nrow(double.idx.mat), 2] < p,
                    double.idx.mat[nrow(double.idx.mat), 1],
                    double.idx.mat[nrow(double.idx.mat), 1] + 1)
    double.idx.mat <- rbind(double.idx.mat, c(first, second))
  }

  # Switch columns in double.idx.mat
  double.idx.mat <- cbind(double.idx.mat[, 2], double.idx.mat[, 1])

  # Compute derivative of transformation
  d.mat <- matrix(NA, nrow = nrow(double.idx.mat), ncol = nrow(double.idx.mat))
  for (i in 1:nrow(double.idx.mat)) {
    for (j in 1:nrow(double.idx.mat)) {
      d.mat[i, j] <- dchol2par.elem(double.idx.mat[j, 1], double.idx.mat[j, 2],
                                    double.idx.mat[i, 1], double.idx.mat[i, 2],
                                    par.chol1)
    }
  }

  # Give appropriate row and column names
  rownames(d.mat) <- apply(double.idx.mat,
                           1,
                           function(row) {sprintf("(%s)", paste(row, collapse= ", "))})
  colnames(d.mat) <- apply(double.idx.mat,
                           1,
                           function(row) {sprintf("(%s)", paste(row, collapse= ", "))})

  # Return the matrix of derivatives
  d.mat
}



#' @title Compute the variance of the estimates.
#'
#' @description This function computes the variance of the estimates computed
#' by the 'estimate.cmprsk.R' function.
#'
#' @param parhatc Vector of estimated parameters, computed in the first part of
#' \code{estimate.cmprsk.R}.
#' @param gammaest Vector of estimated parameters in the regression model for
#' the control function.
#' @param data A data frame.
#' @param admin Boolean value indicating whether the data contains
#' administrative censoring.
#' @param conf Boolean value indicating whether the data contains confounding
#' and hence indicating the presence of \code{z} and, possibly, \code{w}.
#' @param inst Variable encoding which approach should be used for dealing with
#' the confounding. \code{inst = "cf"} indicates that the control function
#' approach should be used. \code{inst = "W"} indicates that the instrumental
#' variable should be used 'as is'. \code{inst = "None"} indicates that Z will
#' be treated as an exogenous covariate. Finally, when \code{inst = "oracle"},
#' this function will access the argument \code{realV} and use it as the values
#' for the control function. Default is \code{inst = "cf"}.
#' @param cf The control function used to estimate the second step.
#' @param eoi.indicator.names Vector of names of the censoring indicator columns
#' pertaining to events of interest. Events of interest will be modeled allowing
#' dependence between them, whereas all censoring events (corresponding to
#' indicator columns not listed in \code{eoi.indicator.names}) will be treated
#' as independent of every other event. If \code{eoi.indicator.names == NULL},
#' all events will be modeled dependently.
#' @param Zbin Indicator value indicating whether (\code{Zbin = TRUE}) or not
#' \code{Zbin = FALSE} the endogenous covariate is binary. Default is
#' \code{Zbin = NULL}, corresponding to the case when \code{conf == FALSE}.
#' @param use.chol Boolean value indicating whether the cholesky decomposition
#' was used in estimating the covariance matrix.
#' @param n.trans Number of competing risks in the model (and hence, number of
#' transformation models).
#' @param totparl Total number of covariate effects (including intercepts) in
#' all of the transformation models combined.
#'
#' @import numDeriv stringr stats
#'
#' @return Variance estimates of the provided vector of estimated parameters.
#'

variance.cmprsk <- function(parhatc, gammaest, data, admin, conf, inst, cf,
                            eoi.indicator.names, Zbin, use.chol, n.trans, totparl) {

  ## Precondition checks

  if (!use.chol) {
    stop("Only implemented for use.chol == TRUE.")
  }
  if (!conf) {
    stop("Only implemented when in cases where there is confounding.")
  }

  ## Define useful variables

  # Number of observations
  n <- nrow(data)

  # Indicator variable whether eoi.indicator.names is specified. If necessary,
  # determine the indices of the events of interest.
  indep <- !is.null(eoi.indicator.names)
  if (indep) {
    cens.ind.columns <- colnames(data)[grepl("delta[1-9][[:digit:]]*", colnames(data))]
    eoi.idxs <- which(cens.ind.columns %in% eoi.indicator.names)
  }

  # Intercept
  intercept <- data[, "x0"]

  # Exogenous covariates
  exo.cov <- data[, grepl("x[1-9][[:digit:]]*", colnames(data)),
                  drop = FALSE]

  # Endogenous covariate and instrumental variable
  Z <- data[, "z"]
  W <- data[, "w"]

  # Design matrix without endogenous covariate
  XandIntercept <- as.matrix(cbind(intercept, exo.cov))
  XandW <- cbind(XandIntercept, W)

  #### Compute the Hessian matrix of the full likelihood at parhatG ####

  # Construct full vector of estimates
  parhatG <- c(parhatc, gammaest)

  # Numerical approximation of the Hessian at parhatG
  Hgamma <- hessian(likIFG.cmprsk.Cholesky,
                    parhatG,
                    data = data,
                    eoi.indicator.names = eoi.indicator.names,
                    admin = admin,
                    conf = conf,
                    Zbin = Zbin,
                    inst = inst,
                    method = "Richardson",
                    method.args = list(eps = 1e-4, d = 0.0001,
                                       zer.tol = sqrt(.Machine$double.eps/7e-7),
                                       r = 6, v = 2, show.details = FALSE)
  )


  # Omit the part of the variance pertaining to the parameters of the control
  # function and compute the inverse.
  H <- Hgamma[1:length(parhatc),1:length(parhatc)]
  HI <- ginv(H)

  #### Compute the variance of estimates ####

  # H_gamma
  Vargamma <- Hgamma[1:length(parhatc),
                     (length(parhatc) + 1):(length(parhatc) + length(gammaest))]

  # Compute all products [1*1, 1*X1, 1*X2, ..., 1*W, X1*X1, X1*X2, ..., X1*W,
  #                       X2*X2, ..., W*W]
  prodvec <- XandW[, 1]
  for (i in 1:length(gammaest)) {
    for (j in 2:length(gammaest)) {
      if (i <= j){
        prodvec <- cbind(prodvec, diag(XandW[, i] %*% t(XandW[, j])))
      }
    }
  }

  # M-matrix: second derivative of m(W, Z, gamma).
  if (Zbin) {

    secder <- t(-dlogis(XandW %*% gammaest)) %*% prodvec

    WM <- matrix(0, nrow = length(gammaest), ncol = length(gammaest))
    WM[lower.tri(WM, diag = TRUE)] <- secder
    WM[upper.tri(WM)] <- t(WM)[upper.tri(WM)]

  } else {

    # Negative of the column sums of prodvec.
    sumsecder <- c(rep(0, ncol(prodvec)))
    for (i in 1:length(sumsecder)) {
      sumsecder[i] <- -sum(prodvec[, i])
    }

    # In this case m(W, Z, gamma) = (Z - W^\top \gamma)². Note that a factor
    # '-2' is omitted.
    WM <- matrix(0, nrow = length(gammaest), ncol = length(gammaest))
    WM[lower.tri(WM, diag = TRUE)] <- sumsecder
    WM[upper.tri(WM)] <- t(WM)[upper.tri(WM)]
  }

  # Inverse of M-matrix
  WMI <- ginv(WM)

  # h_m: first derivative of m(W,Z,gamma). Note that a factor '-2' is omitted.
  mi <- c()
  if (Zbin) {
    diffvec <- Z - plogis(XandW %*% gammaest)
    for(i in 1:n) {
      newrow <- diffvec[i] %*% XandW[i, ]
      mi <- rbind(mi, newrow)
    }
  } else {
    for(i in 1:n) {
      newrow <- cf[i] %*% XandW[i, ]
      mi <- rbind(mi, newrow)
    }
  }
  mi <- t(mi)

  # psi_i-matrix (omitted factors cancel out)
  psii <- -WMI %*% mi

  # h_l(S_i, gamma, delta)
  gi <- c()
  if (indep) {
    for (i in 1:n) {
      J1 <- jacobian(LikI.cmprsk.Cholesky,
                     parhatc,
                     data = data[i, ],
                     eoi.indicator.names = eoi.indicator.names,
                     admin = admin,
                     conf = conf,
                     cf = cf[i],
                     method = "Richardson",
                     method.args = list(eps = 1e-4, d = 0.0001,
                                        zer.tol = sqrt(.Machine$double.eps/7e-7),
                                        r = 6, v = 2, show.details = FALSE))
      gi <- rbind(gi,c(J1))
    }
  } else {
    for (i in 1:n) {
      J1 <- jacobian(likF.cmprsk.Cholesky,
                     parhatc,
                     data = data[i, ],
                     admin = admin,
                     conf = conf,
                     cf = cf[i],
                     method = "Richardson",
                     method.args = list(eps = 1e-4, d = 0.0001,
                                        zer.tol = sqrt(.Machine$double.eps/7e-7),
                                        r = 6, v = 2, show.details = FALSE))
      gi <- rbind(gi,c(J1))
    }
  }

  gi <- t(gi)

  # h_l(S, gamma, delta) + H_gamma %*% Psi_i
  partvar <- gi + Vargamma %*% psii
  Epartvar2 <- (partvar %*% t(partvar))
  totvarex <- HI %*% Epartvar2 %*% t(HI)
  se <- sqrt(abs(diag(totvarex)))

  #### Construct confidence intervals ####

  ## Covariate effect parameters

  # Extract (the variances of) the coefficients. Give appropriate names.
  coeffs <- parhatc[1:totparl]
  coeff.vars <- diag(totvarex[1:totparl, 1:totparl])

  coeff.names <- c()
  for (i in 1:n.trans) {
    for (j in 0:ncol(exo.cov)) {
      coeff.names <- c(coeff.names, sprintf("beta_{%d, %d}", i, j))
    }
    if (conf) {
      coeff.names <- c(coeff.names, sprintf("alpha_{%d}", i))
    }
    if (conf & (inst != "none")) {
      coeff.names <- c(coeff.names, sprintf("lambda_{%d}", i))
    }
  }
  names(coeffs) <- coeff.names
  names(coeff.vars) <- coeff.names

  # Construct the confidence intervals
  coeffs.lb <- coeffs - 1.96 * coeff.vars
  coeffs.ub <- coeffs + 1.96 * coeff.vars

  ## Covariance matrix parameters

  if (indep) {
    s1 <- length(eoi.indicator.names)
  } else {
    s1 <- n.trans
  }

  # Extract Cholesky decomposition parameters
  par.chol <- parhatc[(totparl + 1):(length(parhatc) - n.trans)]
  par.chol1 <- par.chol[1:(choose(s1, 2) + s1)]
  par.chol2 <- NULL
  if ((choose(s1, 2) + s1) < length(par.chol)) {
    par.chol2 <- par.chol[(choose(s1, 2) + s1 + 1):length(par.chol)]
  }

  # Delta method for transforming Cholesky parameter variances to covariance
  # element variances.
  chol.covar.mat <- totvarex[(totparl + 1):(length(parhatc) - n.trans),
                             (totparl + 1):(length(parhatc) - n.trans)]

  chol1.covar.mat <- chol.covar.mat[1:(choose(s1, 2) + s1), 1:(choose(s1, 2) + s1)]
  cov1.elems <- chol2par(par.chol1)
  cov1.vars <- t(dchol2par(par.chol1)) %*% chol1.covar.mat %*% dchol2par(par.chol1)

  if ((choose(s1, 2) + s1) < length(par.chol)) {
    cov2.elems <- par.chol2^2
    chol2.covar.mat <- chol.covar.mat[(choose(s1, 2) + s1 + 1):length(par.chol),
                                      (choose(s1, 2) + s1 + 1):length(par.chol)]
    cov2.vars <- diag(as.matrix(chol2.covar.mat)) * (2 * par.chol2)^2
  }

  # Get all column (row) indices corresponding to variance elements
  var.idxs <- unlist(lapply(names(cov1.elems),
                            function(elem) {
                              var(as.numeric(str_split(gsub("\\(|\\)", "", elem), ", ")[[1]])) == 0
                            }))

  # Extract (the variances of) the variance elements. Give appropriate names.
  var1.elems <- unname(cov1.elems[var.idxs])
  var1.vars <- unname(diag(cov1.vars)[var.idxs])
  var2.elems <- var2.vars <- NULL
  if ((choose(s1, 2) + s1) < length(par.chol)) {
    var2.elems <- cov2.elems
    var2.vars <- cov2.vars
  }

  var.elems <- rep(0, n.trans)
  var.vars <- rep(0, n.trans)
  var.elems[eoi.idxs] <- var1.elems
  var.elems[setdiff(1:n.trans, eoi.idxs)] <- var2.elems
  var.vars[eoi.idxs] <- var1.vars
  var.vars[setdiff(1:n.trans, eoi.idxs)] <- var2.vars

  var.names <- unlist(lapply(1:n.trans, function(i) {sprintf("sigma^2_{%d}", i)}))
  names(var.elems) <- var.names
  names(var.vars) <- var.names

  # Construct the confidence interval for variances in the log-space. Then
  # backtransform.
  log.var.vars <- var.vars / var.elems
  log.var.lb <- log(var.elems) - 1.96 * sqrt(log.var.vars)
  log.var.ub <- log(var.elems) + 1.96 * sqrt(log.var.vars)

  var.lb <- exp(log.var.lb)
  var.ub <- exp(log.var.ub)

  # Extract (the variances of) the covariance elements
  cov1.offdiag.elems <- cov1.elems[!var.idxs]
  cov1.offdiag.vars <- diag(cov1.vars)[!var.idxs]

  # Construct the full covariance matrix (offdiagonal elements)
  cov.mat <- matrix(0, nrow = n.trans, ncol = n.trans)
  cov.mat[eoi.idxs, eoi.idxs][lower.tri(cov.mat[eoi.idxs, eoi.idxs])] <- cov1.offdiag.elems
  cov.mat[upper.tri(cov.mat)] <- t(cov.mat)[upper.tri(cov.mat)]

  # Construct the full variance of covariance matrix (offdiagonal elements)
  cov.var.mat <- matrix(0, nrow = n.trans, ncol = n.trans)
  cov.var.mat[eoi.idxs, eoi.idxs][lower.tri(cov.mat[eoi.idxs, eoi.idxs])] <- cov1.offdiag.vars
  cov.var.mat[upper.tri(cov.var.mat)] <- t(cov.var.mat)[upper.tri(cov.var.mat)]

  # Vectorize the off-diagonal elements. Name the elements
  cov.offdiag.elems <- cov.mat[lower.tri(cov.mat)]
  cov.offdiag.vars <- cov.var.mat[lower.tri(cov.var.mat)]

  cov.offdiag.names <- c()
  for (i in 1:n.trans) {
    for (j in 1:n.trans) {
      if (i < j) {
        cov.offdiag.names <- c(cov.offdiag.names, sprintf("sigma_{%d, %d}", i, j))
      }
    }
  }
  names(cov.offdiag.vars) <- cov.offdiag.names
  names(cov.offdiag.elems) <- cov.offdiag.names

  # Construct the confidence interval for the covariances.
  cov.lb <- cov.offdiag.elems - 1.96 * sqrt(cov.offdiag.vars)
  cov.ub <- cov.offdiag.elems + 1.96 * sqrt(cov.offdiag.vars)

  ## Transformation parameters

  # Extract (the variances of) the transformation parameters. Give appropriate
  # names.
  theta.elems <- parhatc[(length(parhatc) - n.trans + 1):length(parhatc)]
  theta.sds <- se[(length(parhatc) - n.trans + 1):length(parhatc)]

  theta.names <- unlist(lapply(1:n.trans, function(i) {sprintf("theta_{%d}", i)}))
  names(theta.elems) <- theta.names
  names(theta.sds) <- theta.names

  # Confidence interval for transformation parameters
  theta.lb <- theta.elems - 1.96 * theta.sds
  theta.ub <- theta.elems + 1.96 * theta.sds

  # Matrix with all confidence intervals
  CI.mat <- cbind(c(coeffs.lb, var.lb, cov.lb, theta.lb),
                  c(coeffs.ub, var.ub, cov.ub, theta.ub))

  ## Return the results

  # Vector of (transformed) parameters
  params <- c(coeffs, var.elems, cov.offdiag.elems, theta.elems)

  # Vector of variances of (transformed) parameters
  param.vars <- c(coeff.vars, var.vars, cov.offdiag.vars, theta.sds^2)

  # Return results
  list(params = params, param.vars = param.vars, CI.mat = CI.mat)
}



#' @title Print the results from a ptcr model.
#'
#' @description This function takes the output of the function estimate.cmprsk.R
#' as an argument and prints a summary of the results.
#' Note: this function is not yet finished! It should be made a proper S3 method
#' in the future.
#'
#' @param output.estimate.cmprsk The output of the estimate.cmprsk.R function.
#'
#' @noRd
#'
#' @import stringr

mysummary.ptcr <- function(output.estimate.cmprsk) {

  # For now, work with dummy values for the variables
  var.names <- c("(intercept)", "var1", "var2", "var3")
  estimates <- c(1, 20, 300, 4)
  std.error <- c(0.5, 0.141579390248, 200, 4)
  sign <- c("*", "**", ".", "")

  cr.names <- c("cr1", "cr2", "cr3", "cr4")
  cr.vars <- c(1, 0.1234567, 300.1, 4)
  cr.corr <- c(0.1, 0.22, 0.333, 0.4444, 0.55555, 0.666666)
  cr.vars.std <- c(1, 1, 1, 1)
  cr.corr.std <- 1:6

  cr.transfo.param <- c(0, 0.6666666666666, 1.3, 2)

  # Compute the width of each column
  header.names <- c("estimate", "std.error", "sign.")
  col.var.names.width <- max(sapply(var.names, nchar))
  col.estimates.width <- max(sapply(as.character(estimates), nchar), nchar(header.names[1]))
  col.std.error.width <- max(sapply(as.character(std.error), nchar), nchar(header.names[2]))
  col.sign.width <- nchar(header.names[3])

  # Construct header for parameters pertaining to covariates
  header <- paste0(str_pad("", col.var.names.width), "  ",
                   str_pad(header.names[1], col.estimates.width, "right"), "  ",
                   str_pad(header.names[2], col.std.error.width, "right"), "  ",
                   str_pad(header.names[3], col.sign.width, "right"))

  # Construct body for parameters pertaining to covariates
  body <- ""
  for (i in 1:length(estimates)) {
    line <- paste0(str_pad(var.names[i], col.var.names.width, "right"), "  ",
                   str_pad(estimates[i], col.estimates.width, "right"), "  ",
                   str_pad(std.error[i], col.std.error.width, "right"), "  ",
                   str_pad(sign[i], col.sign.width, "right"))

    body <- paste0(body, line)
    if (i < length(estimates)) {
      body <- paste0(body, "\n")
    }
  }

  # Construct header for parameters pertaining to competing risks models

  # Construct body for parameters pertaining to competing risks models

  # Make and print summary
  summary <- paste0(header, "\n", body)
  message(summary)
}




