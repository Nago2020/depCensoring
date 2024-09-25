
#' @title Data generation function for competing risks data
#'
#' @description This function generates competing risk data that can be used in
#' simulation studies.
#'
#' @param n The sample size of the generated data set.
#' @param par List of parameter vectors for each component of the transformation
#' model.
#' @param iseed Random seed.
#' @param s The number of competing risks. Note that the given parameter vector
#' could specify the parameters for the model with more than s competing risks,
#' but in this case only the first s sets of parameters will be considered.
#' @param conf Boolean value indicating whether the data set should contain
#' confounding.
#' @param Zbin Indicator whether the confounded variable is binary
#' \code{Zbin = 1} or not \code{Zbin = 0}. If \code{conf = FALSE}, this
#' variable is ignored.
#' @param Wbin Indicator whether the instrument is binary (\code{Zbin = 1}) or
#' not \code{Zbin = 0}.
#' @param type.cov Vector of characters "c" and "b", indicating which exogenous
#' covariates should be continuous \code{"c"} or binary \code{"b"}.
#' @param A.upper The upper bound on the support of the administrative censoring
#' distribution. This can be used to control for the amount of administrative
#' censoring in the data. Default is \code{A.upper = 15}. \code{A.upper = NULL}
#' will return a data set without administrative censoring.
#'
#' @import mvtnorm MASS
#' @importFrom stats rlogis runif rnorm
#'
#' @return A generated data set
#'

dat.sim.reg.comp.risks = function(n, par, iseed, s, conf, Zbin, Wbin, type.cov,
                                  A.upper = 15) {

  # Check for inconsistencies in the arguments
  n.cov <- length(type.cov) + as.numeric(conf)*2
  n.param <- n.cov + 1
  if (n.param > length(par[[1]])) {
    stop("Too few parameters provided in argument 'par'.")
  }

  # Set randomization seed
  set.seed(iseed)

  # Extract parameters
  sd.idx <- length(par) - 1
  gamma.idx <- length(par)
  beta.mat <- NULL
  for (i in 1:s) {
    beta.mat <- cbind(beta.mat, par[[i]])
  }
  sd <- par[[sd.idx]]
  gamma <- par[[gamma.idx]][1:(nrow(beta.mat) - 1)]

  # If there is no confounding, the parameters for Z and the control function
  # can be ignored.
  if (!conf) {
    beta.mat <- beta.mat[-c(nrow(beta.mat) - 1, nrow(beta.mat)), ]
  }

  # Set some useful parameters
  parl <- nrow(beta.mat)

  # Multivariate normal distribution for error terms
  mu <- rep(0, s)
  sigma <- diag(sd[1:s]^2, nrow = s)
  idx.start <- c(1)
  for (i in 1:(s - 2)) {
    idx.start <- c(idx.start, idx.start[length(idx.start)] + s - i)
  }
  for (i in 1:s) {
    for (j in 1:s) {
      if (i < j) {
        sigma[i, j] <- sd[i]*sd[j]*sd[s + idx.start[i] + j - i - 1]
        sigma[j, i] <- sigma[i, j]
      }
    }
  }
  err <- mvrnorm(n, mu = mu , Sigma = sigma)

  # Intercept
  X <- rep(1, n)

  # Covariates
  for (ct in type.cov) {
    if (ct == "b") {
      X <- cbind(X, sample(c(0, 1), n, replace = TRUE))
    } else {
      X <- cbind(X, rnorm(n = n))
    }
  }
  colnames(X) <- paste0("x", 0:(ncol(X) - 1))

  # If confounded data set, generate confounder and its instrument
  if (conf) {

    # Instrument
    if (Wbin) { # Bernoulli with p = 0.5
      W <- sample(c(0, 1), n, replace = TRUE)
    } else { # Uniform[0, 2]
      W <- runif(n, 0, 2)
    }

    XandW <- as.matrix(cbind(X, W))

    # Confounded variable
    if (Zbin) { # Z is binary
      V <- rlogis(n)
      Z <- as.matrix(as.numeric(XandW %*% gamma - V > 0))
      realV <- (1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
    } else { # Z is continuous
      V <- rnorm(n, 0, 2)
      Z <- XandW %*% gamma + V
      realV <- Z - (XandW %*% gamma)
    }

  }

  # Matrix containing all covariates
  Mgen <- X
  if (conf) {
    Mgen <- cbind(Mgen, Z, realV)
  }

  # Latent times
  times <- matrix(nrow = n, ncol = s)
  for (j in 1:s) {
    times[,j] <- IYJtrans(Mgen %*% beta.mat[,j] + err[,j], sd[length(sd) - s + j])
  }
  if (!is.null(A.upper)) { # administrative censoring
    A <- runif(n, 0, A.upper)
    times <- cbind(times, A)
  }

  # Data matrix
  M <- X
  if (conf) {
    M <- cbind(M, Z, W)
  }

  # Observed time
  Y <- apply(X = times, MARGIN = 1, FUN = min)

  # Censoring indicators
  delta <- matrix(nrow = n, ncol = s)
  for (i in 1:s) {
    delta[,i] <- as.numeric(Y == times[,i])
  }
  if (!is.null(A.upper)) {
    delta <- cbind(delta, as.numeric(Y == times[, s + 1]))
  }

  # data consisting of observed time, censoring indicators, all data and the
  # control function.
  data <- cbind(Y, delta, M)
  if (conf) {
    data <- cbind(data, realV)
  }
  data <- as.data.frame(data)
  cens.names <- paste0("delta", 1:s)
  if (!is.null(A.upper)) {cens.names <- c(cens.names, "da")}
  cov.names <- paste0("X", 0:(ncol(X) - 1))
  if (conf) {conf.names <- c("Z", "W", "realV")} else {conf.names <- NULL}
  colnames(data) <- c("Y", cens.names, cov.names, conf.names)

  # Return the result
  return(data)
}
