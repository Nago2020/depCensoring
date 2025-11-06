

#' @title Obtain vector of dependencies of the QRdepCens code base.
#'
#' @description
#' This function returns a character vector containing all required dependencies
#' for the correct operation of the QRdepCens code base.
#'
#' @noRd
#'
QRdepCens.getDependencies <- function() {
  dep.vec <- c("orthopolynom", "rvinecopulib", "nloptr", "MASS", "numDeriv")
}


#' @title Variance function.
#'
#' @description
#' This function describes the heteroscedasticity in the quantile regression
#' model of D\'Haen et al. (2025).
#'
#' @param X Vector of covariates. Default is \code{X = 1}, representing just
#' an intercept entry (i.e. presenting the homoscedastic case).
#' @param gam Vector of regression coefficients in the model for the variance.
#'
#' @noRd
#'
sigma_fun <- function(X = 1, gam) {
  if (length(X) == 1) {
    return(as.numeric(exp(X %*% gam)))
  } else {
    return(exp(X %*% gam))
  }
}

#' @title Squared L2-norm.
#'
#' @description
#' This function implements the sqaured L2-norm of a vector.
#'
#' @param vect Vector of numerics.
#'
#' @noRd
#'
sqnorm <- function(vect){
  return(sum(vect^2))
}

#' @title Laguerre coefficients
#'
#' @description
#' This function computes the coefficients of the Laguerre polynomial of a
#' Laguerre basis polynomial of a certain degree.
#'
#' @param n Degree of the Laguerre polynomial.
#' @param j Coefficient of \eqn{x^j} in the Laguerre basis polynomial.
#'
#' @noRd
#'
Lag_coeff <- function(n, j){
  choose(n,j)*(-1)^j/factorial(j)
}

#' @title Perturb covariates
#'
#' @description
#' This function perturbs a given vector of covariates using multiplicative
#' uniform noise sampled from \eqn{U[0.5, 2]}. Used to perturb initial values
#' for optimization routine.
#'
#' @param x Vector of covariates.
#'
#' @noRd
#'
perturb <- function(x) {
  return(x*runif(length(x), 0.5, 2))
}

#' @title Slightly perturb covariates
#'
#' @description
#' This function perturbs a given vector of covariates using multiplicative
#' uniform noise sampled from \eqn{U[2/3, 4/3]}. Used to perturb initial values
#' for optimization routine.
#'
#' @param x Vector of numerics.
#'
#' @seealso perturb
#'
#' @noRd
#'
slight_perturb <- function(x) {
  return(x*runif(length(x), 2/3, 4/3))
}

#' @title Partial derivative of Frank copula (wrt second argument).
#'
#' @description
#' This function computes the partial derivative wrt the second argument of a
#' Frank copula, evaluated at a certain point.
#'
#' @param u First argument to partial derivative of copula.
#' @param v Second argument to partial derivative of copula.
#' @param theta Copula parameter.
#'
#' @seealso [h_CT_frank()]
#'
#' @noRd
#'
h_TC_frank <- function(u, v, theta){
  res <- (exp(-theta*v)*(exp(-theta*u) - 1))/
    (exp(-theta) - 1 + (exp(-theta*u) - 1)*(exp(-theta*v) - 1))
  return(pmax(pmin(res, 0.99999),0.00001))
}

#' @title Partial derivative of Frank copula (wrt first argument).
#'
#' @description
#' This function computes the partial derivative wrt the first argument of a
#' Frank copula, evaluated at a certain point.
#'
#' @param u First argument to partial derivative of copula.
#' @param v Second argument to partial derivative of copula.
#' @param theta Copula parameter.
#'
#' @seealso [h_TC_frank()]
#'
#' @noRd
#'
h_CT_frank <- function(u,v,theta){
  res <- (exp(-theta*u)*(exp(-theta*v) - 1))/
    (exp(-theta) - 1 + (exp(-theta*v) - 1)*(exp(-theta*u) - 1))
  return(pmax(pmin(res, 0.99999),0.00001))
}

#' @title Partial derivative of Gumbel copula (wrt second argument).
#'
#' @description
#' This function computes the partial derivative wrt the second argument of a
#' Frank copula, evaluated at a certain point.
#'
#' @param u First argument to partial derivative of copula.
#' @param v Second argument to partial derivative of copula.
#' @param theta Copula parameter.
#'
#' @seealso [h_CT_gumbel()]
#'
#' @noRd
#'
h_TC_gumbel <- function(u,v,theta){
  res <- ((log(u)/log(v))^theta+1)^(-1 + 1/theta)*
    exp(-log(v) -((-log(u))^theta + (-log(v))^theta)^(1/theta))
  return(pmax(pmin(res, 0.99999),0.00001))
}

#' @title Partial derivative of Gumbel copula (wrt first argument).
#'
#' @description
#' This function computes the partial derivative wrt the first argument of a
#' Frank copula, evaluated at a certain point.
#'
#' @param u First argument to partial derivative of copula.
#' @param v Second argument to partial derivative of copula.
#' @param theta Copula parameter.
#'
#' @seealso [h_TC_gumbel()]
#'
#' @noRd
#'
h_CT_gumbel <- function(u,v,theta){
  res <- ((log(v)/log(u))^theta+1)^(-1 + 1/theta)*
    exp(-log(u) -((-log(v))^theta + (-log(u))^theta)^(1/theta))
  return(pmax(pmin(res, 0.99999),0.00001))
}

#' @title Partial derivative of Clayton copula (wrt second argument).
#'
#' @description
#' This function computes the partial derivative wrt the second argument of a
#' Frank copula, evaluated at a certain point.
#'
#' @param u First argument to partial derivative of copula.
#' @param v Second argument to partial derivative of copula.
#' @param theta Copula parameter.
#'
#' @seealso [h_CT_clayton()]
#'
#' @noRd
#'
h_TC_clayton <- function(u,v,theta){
  res <- (1 + (v/u)^theta - v^theta)^(-(theta+1)/theta)
  return(pmax(pmin(res, 0.99999),0.00001))
}

#' @title Partial derivative of Clayton copula (wrt first argument).
#'
#' @description
#' This function computes the partial derivative wrt the first argument of a
#' Frank copula, evaluated at a certain point.
#'
#' @param u First argument to partial derivative of copula.
#' @param v Second argument to partial derivative of copula.
#' @param theta Copula parameter.
#'
#' @seealso [h_TC_clayton()]
#'
#' @noRd
#'
h_CT_clayton <- function(u,v,theta){
  res <- (1 + (u/v)^theta - u^theta)^(-(theta+1)/theta)
  return(pmax(pmin(res, 0.99999),0.00001))
}

#' @title Sample covariates from specified distribution.
#'
#' @description
#' This function samples covariates from a specified distribution. Note that it
#' returns the vector with a prepended intercept term.
#'
#' @param size Number of covariate vectors to sample.
#' @param X_params List of parameters pertaining to the covariate distribution.
#' @param seed Initial seed. Default is \code{seed = NULL}.
#'
#' @noRd
#'
covariate_sampling <- function(size, params_X, seed = NULL) {

  # Set random seed if necessary.
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Extract parameter settings for the covariate distribution.
  X_distribution <- params_X[["X_distribution"]]
  X_dimension <- params_X[["X_dimension"]]
  X_mean <- params_X[["X_mean"]]
  X_var <- params_X[["X_var"]]
  X_lower <- params_X[["X_lower"]]
  X_upper <- params_X[["X_upper"]]

  # Generate the covariates.
  if (X_distribution == "norm") {
    Sigma <- matrix(X_var, nrow = X_dimension)
    return(cbind(1, mvrnorm(n = size, mu = X_mean, Sigma = Sigma)))
  }
  if (X_distribution == "unif") {
    X_dat <- cbind(1, runif(n = size, X_lower, X_upper))
    if (X_dimension == 2) {
      X_dat <- cbind(X_dat, runif(n = size, X_lower, X_upper))
    }
    return(X_dat)
  }
}

#' @title Extract covariates from data frame.
#'
#' @description
#' This function extracts the covariates from the given data frame,
#' corresponding to the column names X1, X2, etc.
#'
#' @param data Data frame.
#'
#' @return Data frame containing the covariates, with column names X1, X2, etc.
#'
#' @noRd
#'
extract.covariates <- function(data) {
  data[, grepl("X[1-9][[:digit:]]*$", colnames(data)), drop = FALSE]
}

#' @title Test if left-hand side is not contained in right-hand side.
#'
#' @description
#' Test whether the elements of the left-hand side vector are contained in the
#' right-hand side vector. This function can be seen as an analogue to the
#' built-in function \code{\%in\%}.
#'
#' @param lhs Left-hand side.
#' @param rhs Right-hand side.
#'
#' @noRd
#'
`%notin%` <- function(lhs, rhs) {
  unlist(lapply(lhs, function(num) {!(num %in% rhs)}))
}

#' @title Check if numeric is integer.
#' @description
#' This function checks if a given numeric value is an integer.
#'
#' @param nbr Numeric.
#'
#' @noRd
#'
is_integer <- function(nbr) {
  nbr == round(nbr)
}

#' @title Compute confidence intervals.
#'
#' @description
#' This function computes 95% confidence intervals based on the provided point
#' estimates and standard deviation estimates.
#'
#' @param est Point estimates.
#' @param sd Standard deviations.
#' @param as.string Boolean flag indicating whether the result should be
#' returned as a sting. Default is \code{as.string = TRUE}.
#'
#' @noRd
#'
make_ci <- function(est, sd, as.string = TRUE) {

  # Construct matrix of confidence intervals.
  ci.mat <- cbind(est - qnorm(0.975) * sd, est + qnorm(0.975) * sd)

  # Construct confidence intervals as strings.
  ci.mat.str <- apply(ci.mat, 1, function(row) {sprintf("[%.2f, %.2f]", round(row[1], 3), round(row[2], 3))})

  # Return the results
  if (as.string) {
    return(ci.mat.str)
  } else {
    return(ci.mat)
  }
}


#' @title Compute initial values for mll optimization for the heteroscedastic
#' model.
#'
#' @description
#' This function determines appropriate starting values for the mll optimization
#' function in the case heteroscedasticity is assumed.
#'
#' @param T1 Subvector of observed times corresponding to \eqn{Delta = 1}, i.e.
#' vector of observed event times.
#' @param X1 Subvector of covariates corresponding to \eqn{Delta = 1}.
#' @param C0 Subvector of observed times corresponding to \eqn{Delta = 0}, i.e.
#' vector of observed censoring times.
#' @param X0 Subvector of covariates corresponding to \eqn{Delta = 1}.
#' @param Delta Vector of censoring indicators.
#' @param cop_name Name of the copula to be used.
#' @param hp List of hyperparameters.
#'
#' @seealso perform_mll_optimisation, determine_init_basis
#'
#' @noRd
#'
determine_init_basis_het <- function(T1, X1, C0, X0, Delta, cop_name, hp) {

  # Extract hyperparameters
  init_value_number_basis <- hp[["init_value_number_basis"]]

  ## data-driven part: beta, gamma and paras_C ##

  # initial parameters for beta and paras_C
  # X already has intercept column: ~ X - 1
  alpha_init <- as.vector(lm(C0 ~ X0 - 1)$coefficients)
  beta_init <- as.vector(lm(T1 ~ X1 - 1)$coefficients)

  # Standard deviation of residuals in linear model for C.
  C_res <- C0 - as.matrix(X0) %*% alpha_init
  sigma_C_init <- sd(C_res)

  # to get an approximate value for gamma
  Z <- (T1 - X1 %*% beta_init)^2

  # avoid near-zero values
  Z[Z < 10^(-100)] <- rep(10^(-100), length((Z[Z < 10^(-100)])))
  gamma_init <- as.vector(lm(log(Z) ~ X1 -1)$coefficients/2)
  gamma_init_update <- c(gamma_init[1] - log(8)/2, gamma_init[-1])

  # used further on in loop for adjusting beta_0
  T_res_no_intercept <- T1 - as.matrix(X1[,-1]) %*% beta_init[-1]

  ## other parameters: random values in admissible/some plausible range
  init_values <- vector(mode = "list", length = init_value_number_basis)

  for (j in 1:init_value_number_basis){
    lambda_init <- runif(n = 1, min = 0.05, max = 0.95)
    transf_lambda_init <- lambda_transform(lambda_init)

    # adjust beta_0 based on the current initial value for lambda
    beta_intercept <- as.numeric(quantile(T_res_no_intercept, lambda_init))
    beta_init_update <- c(beta_intercept, beta_init[-1])

    if (cop_name == "frank"){ # ktau between -1 and 1
      eta_init <- ktau_to_eta(ktau = runif(1, -1, 1), cop_name = cop_name)
    }
    if (cop_name == "frankPos") { #ktau between 0 and 1
      eta_init <- ktau_to_eta(ktau = runif(1, 0, 1), cop_name = cop_name)
    }
    if (cop_name == "gumbel"){ # ktau between 0 and 1
      eta_init <- ktau_to_eta(ktau = runif(1,0,1), cop_name = cop_name)
    }
    if (cop_name == "clayton"){ # ktau between 0 and 1
      eta_init <- ktau_to_eta(ktau = runif(1,0,1), cop_name = cop_name)
    }
    if (cop_name == "indep"){
      eta_init <- 0 # not used during optimisation anyway
    }

    init_values[[j]] <- c(eta_init, perturb(beta_init_update),
                          perturb(gamma_init_update),
                          transf_lambda_init,
                          perturb(alpha_init),
                          log(perturb(sigma_C_init)))
  }
  return(init_values)
}

#' @title Compute initial values for mll optimization (homoscedastic case).
#'
#' @description
#' This function determines appropriate starting values for the mll optimization
#' function in the case that homoscedasticity can be assumed.
#'
#' @param T1 Subvector of observed times corresponding to \eqn{Delta = 1}, i.e.
#' vector of observed event times.
#' @param X1 Subvector of covariates corresponding to \eqn{Delta = 1}.
#' @param C0 Subvector of observed times corresponding to \eqn{Delta = 0}, i.e.
#' vector of observed censoring times.
#' @param X0 Subvector of covariates corresponding to \eqn{Delta = 1}.
#' @param Delta Vector of censoring indicators.
#' @param cop_name Name of the copula to be used.
#' @param hp List of hyperparameters.
#'
#' @seealso perform_mll_optimisation, determine_init_basis
#'
#' @noRd
#'
determine_init_basis_hom <- function(T1, X1, C0, X0, Delta, cop_name, hp) {

  # Extract hyperparameters
  init_value_number_basis <- hp[["init_value_number_basis"]]
  lambda <- hp[["lambda"]]

  ## data-driven part: beta, gamma and paras_C ##

  # initial parameters for beta and paras_C
  # X already has intercept column: ~ X - 1
  alpha_init <- as.vector(lm(C0 ~ X0 - 1)$coefficients)
  beta_init <- as.vector(lm(T1 ~ X1 - 1)$coefficients)

  # adjust beta_0 based on current initial values & (fixed) lambda
  T_res_no_intercept <- T1 - as.matrix(X1[,-1]) %*% beta_init[-1]
  beta_intercept <- as.numeric(quantile(T_res_no_intercept, lambda))
  beta_init_update <- c(beta_intercept, beta_init[-1])

  C_res <- C0 - as.matrix(X0) %*% alpha_init
  sigma_C_init <- sd(C_res)


  # initial parameter for gamma_T
  sigma_T_init <- sqrt(var(T_res_no_intercept) *
                         (lambda^2*(1 - lambda)^2)/(lambda^2 + (1-lambda)^2))
  gamma_T_init <- log(sigma_T_init)

  ## remaining parameter ktau/eta: random values in admissible range
  init_values <- vector(mode = "list", length = init_value_number_basis)

  for (j in 1:init_value_number_basis){
    if (cop_name == "frank"){ # ktau between -1 and 1
      eta_init <- ktau_to_eta(ktau = runif(1, -1, 1), cop_name = cop_name)
    }
    if (cop_name == "frankPos") { #ktau between 0 and 1
      eta_init <- ktau_to_eta(ktau = runif(1, 0, 1), cop_name = cop_name)
    }
    if (cop_name == "gumbel"){ # ktau between 0 and 1
      eta_init <- ktau_to_eta(ktau = runif(1,0,1), cop_name = cop_name)
    }
    if (cop_name == "clayton"){ # ktau between 0 and 1
      eta_init <- ktau_to_eta(ktau = runif(1,0,1), cop_name = cop_name)
    }
    if (cop_name == "indep"){
      eta_init <- 0 # not used during optimisation anyway
    }

    init_values[[j]] <- c(eta_init, perturb(beta_init_update),
                          perturb(gamma_T_init),
                          perturb(alpha_init),
                          log(perturb(sigma_C_init)))
  }
  return(init_values)
}

#' @title Compute initial values for mll optimization.
#'
#' @description
#' This function determines appropriate starting values for the mll optimization
#' function. This function is simply a wrapper for the homoscedastic and
#' heteroscedastic case.
#'
#' @param T1 Subvector of observed times corresponding to \eqn{Delta = 1}, i.e.
#' vector of observed event times.
#' @param X1 Subvector of covariates corresponding to \eqn{Delta = 1}.
#' @param C0 Subvector of observed times corresponding to \eqn{Delta = 0}, i.e.
#' vector of observed censoring times.
#' @param X0 Subvector of covariates corresponding to \eqn{Delta = 1}.
#' @param Delta Vector of censoring indicators.
#' @param cop_name Name of the copula to be used.
#' @param hp List of hyperparameters.
#'
#' @seealso perform_mll_optimisation, determine_init_basis_hom,
#' determine_init_basis_het
#'
#' @noRd
#'
determine_init_basis <- function(T1, X1, C0, X0, Delta, cop_name, hp) {
  if (hp$homoscedastic) {
    determine_init_basis_hom(T1, X1, C0, X0, Delta, cop_name, hp)
  } else {
    determine_init_basis_het(T1, X1, C0, X0, Delta, cop_name, hp)
  }
}

#' @title Determine list of initial values.
#'
#' @description
#' This function determines a list of initial values based on slightly
#' perturbing the given vector of parameters.
#'
#' @param current_paras Parameter vector. In the implementation, this is the
#' intermediate result after degree selection of the Laguerre polynomials.
#' @param hp List of hyperparameters.
#' @param indep_assumption Boolean flag indicating whether independence can be
#' assumed. Default is \code{indep_assumption = FALSE}.
#'
#' @noRd
#'
determine_init_final <- function(current_paras, hp, indep_assumption = FALSE) {
  init_value_number_final <- hp[["init_value_number_final"]]
  new_init_values <- vector(mode = "list", length = init_value_number_final)
  if (!indep_assumption) {
    for (j in seq_len(init_value_number_final)) {
      new_init_values[[j]] <- slight_perturb(x = current_paras)
    }
  } else {
    for (j in seq_len(init_value_number_final)) {
      new_init_values[[j]] <- c(0, slight_perturb(x = current_paras))
    }
  }
  return(new_init_values)
}


#' @title Density function of error term for C (homoscedastic + Gaussian case).
#'
#' @description
#' This function defines the density function of the error term in the linear
#' model for C, in the homoscedastic case and when a normal distribution is
#' assumed on C.
#'
#' @param y Outcome variables; i.e. censoring times.
#' @param x Predictors; i.e. covariates.
#' @param alpha_C Regression coefficients.
#' @param sigma_C Standard deviation.
#'
#' @noRd
#'
dCens_norm <- function(y, x, alpha_C, sigma_C){
  z <- 1/sigma_C*(y -  x %*% alpha_C)
  return(1/sigma_C*dnorm(z))
}

#' @title Distribution function of error term for C (homoscedastic + Gaussian
#' case).
#'
#' @description
#' This function defines the distribution function of the error term in the linear
#' model for C, in the homoscedastic case and when a normal distribution is
#' assumed on C.
#'
#' @inheritParams dCens_norm
#'
#' @noRd
#'
pCens_norm <- function(y, x, alpha_C, sigma_C){
  z <- 1/sigma_C*(y -  x %*% alpha_C)
  return(pnorm(z))
}

#' @title Quantile function of error term for C (homoscedastic + Gaussian case).
#'
#' @description
#' This function defines the quantile function of the error term in the linear
#' model for C, in the homoscedastic case and when a normal distribution is
#' assumed on C.
#'
#' @inheritParams dCens_norm
#'
#' @noRd
#'
qCens_norm <- function(w, x, alpha_C, sigma_C) {
  return(qnorm(w)*sigma_C + x %*% alpha_C)
}

#' @title Density function of error term for C (homoscedastic case).
#'
#' @description
#' This function defines the density function of the error term in the linear
#' model for C, in the homoscedastic case.
#'
#' @param y Outcome variables; i.e. censoring times.
#' @param x Predictors; i.e. covariates.
#' @param distr_C Name of distribution for C.
#' @param params_C List of parameters for the model for C. Default is
#' \code{params_C = NULL}. However, a \code{NULL} value is only valid in the
#' specific case where the function is called with only three arguments.
#'
#' @noRd
#'
dCens <- function(y, x, distr_C, params_C = NULL) {

  # -- For compatibility with estimation methodology --
  # Throughout the implementation of the estimation procedure, this function
  # is used with arguments dCens(y, x, transf_params_C) and the distribution of
  # C is assumed to be normal.
  if (is.null(params_C) & !is.character(distr_C)) {
    transf_paras_C <- distr_C
    alpha_C <- transf_paras_C[-length(transf_paras_C)]
    sigma_C <- exp(transf_paras_C[length(transf_paras_C)])
    distr_C <- "norm"
    params_C <- list(alpha_C = alpha_C, sigma_C = sigma_C)
  }

  # If C follows a normal distribution...
  if (distr_C == "norm") {
    alpha_C <- params_C[["alpha_C"]]
    sigma_C <- params_C[["sigma_C"]]
    return(dCens_norm(y, x, alpha_C, sigma_C))
  }

  # If non of the above distribution applied, throw an error.
  stop("Distribution for C not implemented yet.")
}

#' @title Distribution function of error term for C (homoscedastic case).
#'
#' @description
#' This function defines the distribution function of the error term in the linear
#' model for C, in the homoscedastic case.
#'
#' @inheritParams dCens
#'
#' @noRd
#'
pCens <- function(y, x, distr_C, params_C = NULL) {

  # -- For compatibility with estimation methodology --
  # Throughout the implementation of the estimation procedure, the this function
  # is used with arguments dCens(y, x, transf_params_C) and the distribution of
  # C is assumed to be normal.
  if (is.null(params_C) & !is.character(distr_C)) {
    transf_paras_C <- distr_C
    alpha_C <- transf_paras_C[-length(transf_paras_C)]
    sigma_C <- exp(transf_paras_C[length(transf_paras_C)])
    distr_C <- "norm"
    params_C <- list(alpha_C = alpha_C, sigma_C = sigma_C)
  }

  # If C follows a normal distribution...
  if (distr_C == "norm") {
    alpha_C <- params_C[["alpha_C"]]
    sigma_C <- params_C[["sigma_C"]]
    return(pCens_norm(y, x, alpha_C, sigma_C))
  }

  # If non of the above distribution applied, throw an error.
  stop("Distribution for C not implemented yet.")
}

#' @title Quantile function of error term for C (homoscedastic case).
#'
#' @description
#' This function defines the quantile function of the error term in the linear
#' model for C, in the homoscedastic case.
#'
#' @inheritParams dCens
#'
#' @noRd
#'
qCens <- function(y, x, distr_C, params_C = NULL) {

  # -- For compatibility with estimation methodology --
  # Throughout the implementation of the estimation procedure, the this function
  # is used with arguments dCens(y, x, transf_params_C) and the distribution of
  # C is assumed to be normal.
  if (is.null(params_C) & !is.character(distr_C)) {
    transf_paras_C <- distr_C
    alpha_C <- transf_paras_C[-length(transf_paras_C)]
    sigma_C <- exp(transf_paras_C[length(transf_paras_C)])
    distr_C <- "norm"
    params_C <- list(alpha_C = alpha_C, sigma_C = sigma_C)
  }

  # If C follows a normal distribution...
  if (distr_C == "norm") {
    alpha_C <- params_C[["alpha_C"]]
    sigma_C <- params_C[["sigma_C"]]
    return(qCens_norm(y, x, alpha_C, sigma_C))
  }

  # If non of the above distribution applied, throw an error.
  stop("Distribution for C not implemented yet.")
}

#### 4)  Transformation functions ####

#' @title Compute Kendall's tau based on \eqn{eta}.
#'
#' @param eta Transformed Kendall's tau parameter, which can take values on the
#' entire real line.
#' @param cop_name Name of copula under consideration. Could be either
#' \code{"frank"}, \code{"frankPos"}, \code{"gumbel"}, \code{"clayton"} or
#' \code{"indep"}.
#'
#' @returns Kendall's tau.
#'
#' @noRd
#'
eta_to_ktau <- function(eta, cop_name){
  if (cop_name == "frank"){ # ktau between -1 and 1
    ktau <- tanh(eta)
  }
  if (cop_name %in% c("frankPos", "gumbel", "clayton")) { #ktau between 0 and 1
    ktau <- 1/(1 + exp(-eta))
  }
  if (cop_name == "indep"){ # no ktau present
    ktau <- 0
  }
  return(ktau)
}

#' @title Compute \eqn{eta} based on Kendall's tau.
#'
#' @param ktau Value of Kendall's tau.
#' @param cop_name Name of copula under consideration. Could be either
#' \code{"frank"}, \code{"frankPos"}, \code{"gumbel"}, \code{"clayton"} or
#' \code{"indep"}.
#'
#' @returns eta.
#'
#' @noRd
#'
ktau_to_eta <- function(ktau, cop_name){
  if (cop_name == "frank"){ # ktau between -1 and 1
    eta <- 0.5* log((1 + ktau)/(1 - ktau))
  }
  if (cop_name %in% c("frankPos", "gumbel", "clayton")) { # ktau between 0 and 1
    eta <- log(ktau/(1 - ktau))
  }
  if (cop_name == "indep"){ # no ktau present
    eta <- 0
  }
  return(eta)
}

#' @title Compute copula parameter based on Kendall's tau.
#'
#' @param ktau Value of Kendall's tau.
#' @param cop_name Name of copula under consideration. Could be either
#' \code{"frank"}, \code{"frankPos"}, \code{"gumbel"}, \code{"clayton"} or
#' \code{"indep"}.
#'
#' @importFrom rvinecopulib ktau_to_par
#'
#' @returns Copula parameter.
#'
#' @noRd
#'
ktau_to_theta <- function(ktau, cop_name) {
  cop_name_rvcl <- cop_name
  if (cop_name == "frankPos") {
    cop_name_rvcl <- "frank"
  }
  theta <- as.numeric(rvinecopulib::ktau_to_par(family = cop_name_rvcl, tau = ktau))
  return(theta)
}

#' @title Transform quantile value to whole real line.
#'
#' @description
#' This function transforms the quantile value, which is probability and hence
#' a value between [0, 1] to the whole real line.
#'
#' @param lambda Quantile value.
#'
#' @returns Transformed quantile value.
#'
#' @seealso lambda_inverse_transform.
#'
#' @noRd
#'
lambda_transform <- function(lambda){
  return(log(lambda/(1 - lambda)))
}

#' @title Back-transform quantile value
#'
#' @description
#' This function transforms the quantile value represented as a number on the
#' whole real line to one that takes values in [0, 1].
#'
#' @param transf_lambda Transformed quantile value.
#'
#' @returns Back-transformed quantile value.
#'
#' @seealso lambda_transform.
#'
#' @noRd
#'
lambda_inverse_transform <- function(transf_lambda){
  return(1/(1 + exp(-transf_lambda)))
}

#' @title Back-transform the parameter vector.
#'
#' @description
#' This function back-transforms certain elements of the parameter vector. In
#' particular, it targets the dependence parameter, variance parameter and - if
#' applicable - lambda. I.e., all parameters that are restricted to some proper
#' subset of the real line.
#'
#' @param transf_para_vec Transformed parameter vector.
#' @param cop_name Name of considered copula.
#' @param hp List of hyperparameters.
#'
#' @returns Transformed parameter vector.
#'
#' @noRd
#'
para_inv_transform <- function(transf_para_vec, cop_name, hp) {

  # Extract hyperparameter
  lambda_index <- hp[["lambda_index"]]
  homoscedastic <- hp[["homoscedastic"]]

  para_vec <- transf_para_vec

  if (cop_name == "indep"){
    ktau_length <- 0
  } else {
    ktau_length <- 1
  }

  # replace eta by ktau, if present
  if (ktau_length == 1){
    para_vec[1] <- eta_to_ktau(eta = para_vec[1], cop_name = cop_name)
  }

  # replace transf_lambda by lambda
  if (!homoscedastic) {
    para_vec[lambda_index + ktau_length] <- lambda_inverse_transform(
      para_vec[lambda_index + ktau_length])
  }

  # replace log_sigma_C by sigma_C
  para_vec[length(para_vec)] <- exp(para_vec[length(para_vec)])

  return(para_vec)
}


#' @title Density function of EAL distribution.
#'
#' @description
#' This function implements the density function of an Enriched Asymetric
#' Laplace (EAL) distribution.
#'
#' @param y (vector of) Observed time(s).
#' @param x Covariate vector/matrix.
#' @param beta Coefficients in linear quantile regression model.
#' @param gam_par Parameters for variance model (heteroscedastic case.)
#' @param lambda Quantile of interest.
#' @param phi_til Laguerre expansion coefficients for y <= 0.
#' @param phi Laguerre expansion coefficients for y > 0.
#' @param hp List of hyperparameters.
#'
#' @returns EAL density function evaluation.
#'
#' @noRd
#'
dEAL <- function(y, x, beta, gam_par, lambda, phi_til, phi, hp) {

  # Standard deviation in quantile regression model for T.
  x.simga_fun <- if (hp$homoscedastic) {1} else {x}
  s <- sigma_fun(x.simga_fun, gam_par)

  # Point at which to evaluate the EAL distribution; residuals in linear
  # quantile regression model.
  z <- 1/s * (y - x %*% beta)

  # Laguerre polynomial coefficients.
  phi_til_ext <- c(1, phi_til)
  phi_ext <- c(1, phi)

  # Lag of degree 0 is constant 1, so first term in sum is always fixed constant
  left_sum <- 1
  right_sum <- 1

  poly_list <- orthopolynom::laguerre.polynomials(max(length(phi),length(phi_til)),
                                                  normalized = TRUE)

  if (length(y) == 1){ # avoid indexing one vector x as x[pos_indices,]
    if (z <= 0){
      for (k in seq_along(phi_til)){
        poly <- as.function(poly_list[[k+1]])
        left_sum <- left_sum + phi_til[k]*poly((lambda-1)*z)
      }
      result <- (1-lambda)*lambda/s*1/sqnorm(phi_til_ext)*exp(-z *(lambda-1)) * left_sum^2

    } else {
      for (k in seq_along(phi)){
        poly <- as.function(poly_list[[k+1]])
        right_sum <- right_sum + phi[k]*poly(lambda*z)
      }
      result <- lambda*(1-lambda)/s*1/sqnorm(phi_ext)*exp(-z*lambda) * right_sum^2
    }
  } else {
    pos_indices <- which(z > 0)

    z_pos <- z[pos_indices]
    z_neg <- z[-pos_indices]

    if (hp$homoscedastic) {
      s_pos <- s
      s_neg <- s
    } else {
      s_pos <- s[pos_indices]
      s_neg <- s[-pos_indices]
    }

    for (k in seq_along(phi_til)){
      poly <- as.function(poly_list[[k+1]])
      left_sum <- left_sum + phi_til[k]*poly((lambda-1)*z_neg)
    }

    for (k in seq_along(phi)){
      poly <- as.function(poly_list[[k+1]])
      right_sum <- right_sum + phi[k]*poly(lambda*z_pos)
    }

    result <- rep(0, length(z))
    result[pos_indices] <- lambda*(1-lambda)/s_pos*1/sqnorm(phi_ext)*
      exp(-z_pos*lambda) * right_sum^2
    result[-pos_indices] <- (1-lambda)*lambda/s_neg*1/sqnorm(phi_til_ext)*
      exp(-z_neg *(lambda-1)) * left_sum^2
  }
  return(as.matrix(result))
}

#' @title Left side of EAL distribution function.
#'
#' @description
#' This function implements the left side of the distribution function of an
#' Enriched Asymetric Laplace (EAL) distribution.
#'
#' @param y (vector of) Observed time(s).
#' @param x Covariate vector/matrix.
#' @param beta Coefficients in linear quantile regression model.
#' @param gam_par Parameters for variance model (heteroscedastic case.)
#' @param lambda Quantile of interest.
#' @param phi_til Laguerre expansion coefficients for y <= 0.
#' @param phi Laguerre expansion coefficients for y > 0.
#'
#' @returns EAL distribution function evaluation.
#'
#' @seealso pEAL
#'
#' @noRd
#'
pEAL_neg <- function(y, x, beta, gam_par, lambda, phi_til) {

  # Determine if homoscedasticity is assumed
  homoscedastic <- length(gam_par) == 1

  # Standard deviation in quantile regression model for T.
  x.simga_fun <- if (homoscedastic) {1} else {x}
  s <- sigma_fun(x.simga_fun, gam_par)

  z <- 1/s * (y - x %*% beta)

  phi_til_ext <- c(1, phi_til)
  n <- length(phi_til)

  sum <- 0
  for (a in 0:n){
    for (b in 0:n){
      for (c in 0:a){
        for (d in 0:b){
          inner_sum <- 0
          for (j in 0:(c+d)){
            inner_sum <- inner_sum + factorial(c+d)/factorial(j)*
              (-(1-lambda)*z)^j
          }
          sum <- sum + phi_til_ext[a+1]*phi_til_ext[b+1] *
            Lag_coeff(a,c)*Lag_coeff(b,d)* inner_sum
        }
      }
    }
  }
  res <- lambda*exp((1-lambda)*z) * sum/sqnorm(phi_til_ext)

  res_corrected <- pmin(pmax(0, res), 1)
  return(res_corrected) # could still be NaN !
}

#' @title Right side of EAL distribution function.
#'
#' @description
#' This function implements the right side of the distribution function of an
#' Enriched Asymetric Laplace (EAL) distribution.
#'
#' @param y (vector of) Observed time(s).
#' @param x Covariate vector/matrix.
#' @param beta Coefficients in linear quantile regression model.
#' @param gam_par Parameters for variance model (heteroscedastic case.)
#' @param lambda Quantile of interest.
#' @param phi_til Laguerre expansion coefficients for y <= 0.
#' @param phi Laguerre expansion coefficients for y > 0.
#'
#' @returns EAL distribution function evaluation.
#'
#' @seealso pEAL
#'
#' @noRd
#'
pEAL_pos <- function(y, x, beta, gam_par, lambda, phi) {

  # Determine if homoscedasticity is assumed
  homoscedastic <- length(gam_par) == 1

  # Standard deviation in quantile regression model for T.
  x.simga_fun <- if (homoscedastic) {1} else {x}
  s <- sigma_fun(x.simga_fun, gam_par)

  z <- 1/s * (y - x %*% beta)

  phi_ext <- c(1, phi)
  n <- length(phi)

  sum <- 0
  for (a in 0:n){
    for (b in 0:n){
      for (c in 0:a){
        for (d in 0:b){
          inner_sum <- 0
          for (j in 0:(c+d)){
            inner_sum <- inner_sum + factorial(c+d)/factorial(j)*(lambda*z)^j
          }
          sum <- sum + phi_ext[a+1]*phi_ext[b+1]*Lag_coeff(a,c)*Lag_coeff(b,d)*
            (factorial(c+d)- exp(-lambda*z)*inner_sum)
        }
      }
    }
  }
  res <- lambda + (1-lambda) * sum/sqnorm(phi_ext)

  res_corrected <- pmin(pmax(0, res), 1)
  return(res_corrected) # could still be NaN !
}

#' @title EAL distribution function.
#'
#' @description
#' This function implements the distribution function of an Enriched Asymetric
#' Laplace (EAL) distribution. It is mainly a wrapper around the left and right
#' side distribution function.
#'
#' @param y (vector of) Observed time(s).
#' @param x Covariate vector/matrix.
#' @param beta Coefficients in linear quantile regression model.
#' @param gam_par Parameters for variance model (heteroscedastic case.)
#' @param lambda Quantile of interest.
#' @param phi_til Laguerre expansion coefficients for y <= 0.
#' @param phi Laguerre expansion coefficients for y > 0.
#'
#' @returns EAL distribution function evaluation.
#'
#' @seealso [pEAL_neg(), pEAL_pos()]
#'
#' @noRd
#'
pEAL <- function(y, x, beta, gam_par, lambda, phi_til, phi) {

  # Determine if homoscedasticity is assumed
  homoscedastic <- length(gam_par) == 1

  # Standard deviation in quantile regression model for T.
  x.simga_fun <- if (homoscedastic) {1} else {x}
  s <- sigma_fun(x.simga_fun, gam_par)

  z <- 1/s * (y - x %*% beta)


  if (length(y) == 1){ # avoid indexing one vector x as x[pos_indices,]
    if (z < 0){
      res_vec <- pEAL_neg(y, x, beta, gam_par, lambda, phi_til)
    } else {
      res_vec <- pEAL_pos(y, x, beta, gam_par, lambda, phi)
    }
  } else { # general case, correct subsetting dimensions
    pos_indices <- which(z > 0)
    res_vec <- as.matrix(rep(0, length(z)))

    if (length(pos_indices) == 0){
      res_vec <- pEAL_neg(y, x, beta, gam_par, lambda, phi_til)
    } else {
      res_vec[pos_indices,] <- pEAL_pos(y[pos_indices], x[pos_indices,],
                                        beta, gam_par, lambda, phi)


      if (length(pos_indices) < length(z)){
        res_vec[-pos_indices,] <- pEAL_neg(y[-pos_indices], x[-pos_indices,],
                                           beta, gam_par, lambda, phi_til)
      }
    }
  }

  # correct for numerical errors that lead to values outside [0,1]
  res_vec <- pmin(pmax(0, res_vec), 1)

  return(as.matrix(res_vec)) # could still be NaN !
}

#' @title EAL quantile function.
#'
#' @description
#' This function implements the quantile function of an Enriched Asymetric
#' Laplace (EAL) distribution.
#'
#' @param w (vector of) Observed probability/ies.
#' @param x Covariate vector/matrix.
#' @param beta Coefficients in linear quantile regression model.
#' @param gam_par Parameters for variance model (heteroscedastic case.)
#' @param lambda Quantile of interest.
#' @param phi_til Laguerre expansion coefficients for y <= 0.
#' @param phi Laguerre expansion coefficients for y > 0.
#'
#' @returns EAL quantile function evaluation.
#'
#' @noRd
#'
qEAL <- function(w, x, beta, gam_par, lambda, phi_til, phi) {

  # Set some useful parameters
  X_dimension <- length(beta) - 1

  # Determine if homoscedasticity is assumed
  homoscedastic <- length(gam_par) == 1

  # determine basis quantile from EAL, only afterwards shift
  x_null = matrix(rep(1, length(w)*(X_dimension + 1)), nrow = length(w))
  beta_null = rep(0, X_dimension + 1)
  gam_par_null = if (homoscedastic) {0} else {rep(0, X_dimension + 1)}

  f_null <- function(y){
    return(pEAL(y, x_null, beta_null, gam_par_null, lambda, phi_til, phi) - w)
  }

  basis_quantile <- uniroot(f_null, lower = -10, upper = 10, extendInt = "yes")$root

  x.sigma_fun <- if (homoscedastic) {1} else {x}
  s <- sigma_fun(x.sigma_fun, gam_par)

  return(x %*% beta + s*basis_quantile)
}



#' @title Log-likelihood function of D'Haen et al. (2025).
#'
#' @description This function implements the log-likelihood function used in
#' the method of D'Haen et al. (2025). In particular, it implements equation
#' (4.1) of their paper.
#'
#' @param all_paras Vector containing all parameters.
#' @param lag_degs Degrees of Laguerre polynomials.
#' @param cop_name Name of considered copula.
#' @param Y Vector of observed times.
#' @param X Design matrix.
#' @param Delta Vector of censoring indicators.
#' @param hp List of hyperparameters.
#'
#' @returns Log-likelihood evaluation.
#'
#' @noRd
#'
#' @references D'Haen, M., Van Keilegom, I. and Verhasselt, A. (2025). Quantile
#' regression under dependent censoring with unknown association. Lifetime Data
#' Analysis 31(2):253-299.
#'
log_lik_fun <- function(all_paras, lag_degs, cop_name, Y, X, Delta, hp) {

  # 0) Extract hyperparameters

  beta_full_index <- hp[["beta_full_index"]]
  gamma_index <- hp[["gamma_index"]]
  lambda_index <- hp[["lambda_index"]]
  phi_til_index <- hp[["phi_til_index"]]
  phi_index <- hp[["phi_index"]]
  lambda <- hp[["lambda"]]
  homoscedastic <- hp[["homoscedastic"]]
  heteroscedastic <- !homoscedastic

  # 1) Retrieving current parameters & corresponding copula object

  paras_T_indices <- c(beta_full_index, gamma_index, lambda_index,
                       phi_til_index(lag_degs),
                       phi_index(lag_degs)) + 1

  eta <- all_paras[1]
  transf_paras_T <- all_paras[paras_T_indices]
  transf_paras_C <- all_paras[-c(1, paras_T_indices)]

  beta_full <- transf_paras_T[beta_full_index]
  if (heteroscedastic) {
    transf_lam <- transf_paras_T[lambda_index]
    lambda <- lambda_inverse_transform(transf_lam)
  }
  gam_par <- transf_paras_T[gamma_index]
  phi_til <- transf_paras_T[phi_til_index(lag_degs)]
  phi <- transf_paras_T[phi_index(lag_degs)]

  # 2) Construct argument matrices for censored & uncensored cases
  Y_un <- Y[Delta == 1]
  Y_cens <- Y[Delta == 0]

  X_un <- X[Delta == 1, ]
  X_cens <- X[Delta == 0, ]

  # 3) Density, distribution & copula partial derivative evaluation

  ## density
  dens_T_un <- dEAL(Y_un, X_un, beta_full, gam_par, lambda, phi_til, phi, hp)
  dens_C_cens <- dCens(Y_cens, X_cens, transf_paras_C)

  ## distribution
  distr_T_un <- pEAL(Y_un, X_un, beta_full, gam_par, lambda, phi_til, phi)
  distr_T_cens <- pEAL(Y_cens, X_cens, beta_full, gam_par, lambda, phi_til, phi)

  distr_C_un <- pCens(Y_un, X_un, transf_paras_C)
  distr_C_cens <- pCens(Y_cens, X_cens, transf_paras_C)


  ## copula partial derivatives
  theta <- ktau_to_theta(eta_to_ktau(eta, cop_name), cop_name)

  if (cop_name %in% c("frank", "frankPos")) {
    h_C_given_T_uncens <- h_CT_frank(u = distr_T_un, v = distr_C_un, theta)
    h_T_given_C_cens <- h_TC_frank(u = distr_T_cens, v = distr_C_cens, theta)
  }
  if (cop_name == "gumbel"){
    h_C_given_T_uncens <- h_CT_gumbel(u = distr_T_un, v = distr_C_un, theta)
    h_T_given_C_cens <- h_TC_gumbel(u = distr_T_cens, v = distr_C_cens, theta)
  }
  if (cop_name == "clayton"){
    h_C_given_T_uncens <- h_CT_clayton(u = distr_T_un, v = distr_C_un, theta)
    h_T_given_C_cens <- h_TC_clayton(u = distr_T_cens, v = distr_C_cens, theta)
  }

  # Return the results

  lh_uncens <- dens_T_un * (1 - h_C_given_T_uncens)
  lh_cens <- dens_C_cens * (1 - h_T_given_C_cens)

  llh <- sum(log(lh_uncens)) + sum(log(lh_cens))

  res <- min(max(llh, -100000), 100000) # avoid returning +- Inf
  if (is.nan(res)){ # avoid NaN (due to Inf - Inf by num. errors) as well
    res <- -100000
  }
  return(res)
}

#' @title Negative log-likelihood function of D'Haen et al. (2025).
#'
#' @description This function implements the log-likelihood function used in
#' the method of D'Haen et al. (2025). In particular, it implements equation
#' (4.1) of their paper. It will return the negative of the
#' log-likelihood. This function is simply a wrapper around
#' \code{log_lik_fun}.
#'
#' @inheritParams log_lik_fun
#' @returns Negative of the log-likelihood evaluation
#' @seealso log_lik_fun
#'
#' @noRd
#'
neg_log_lik_fun <- function(all_paras, lag_degs, cop_name, Y, X, Delta, hp) {
  -log_lik_fun(all_paras, lag_degs, cop_name, Y, X, Delta, hp)
}

#' @title Continuity constraint of EAL distribution.
#'
#' @description
#' This function implements the equality constraint (6.1) of D'Haen et al.
#' (2025) that ensures continuity of the EAL distribution. It is implemented
#' as a combination of two inequality constraints.
#'
#' @param all_paras Vector containing all parameters.
#' @param lag_degs Laguerre polynomial maximal degrees.
#' @param cop_name Name of considered copula.
#' @param Y Vector of observed times.
#' @param X Design matrix.
#' @param Delta Vector of censoring indicators.
#' @param hp List of hyperparameters.
#'
#' @returns Vector of continuity constraint violations.
#'
#' @noRd
#'
ctuity_constr_fun <- function(all_paras, lag_degs,
                              cop_name, Y, X, Delta, hp) {

  # Extract hyperparameters
  phi_til_index <- hp[["phi_til_index"]]
  phi_index <- hp[["phi_index"]]

  # args. on second line not needed, but should be the same as for obj. fun.

  phi_til <- all_paras[phi_til_index(lag_degs) + 1] # omit eta
  phi <- all_paras[phi_index(lag_degs) + 1] # omit eta

  phi_til_ext <- c(1, phi_til)
  phi_ext <- c(1, phi)

  # combination of Lag. polynomials evaluated in 0 on neg./pos. axis:
  left <- sum(phi_til_ext)^2 / sqnorm(phi_til_ext)
  right <- sum(phi_ext)^2 / sqnorm(phi_ext)

  jump <- left - right

  return(c(jump, - jump))
}

#' @title Optimize the log-likelihood function.
#'
#' @description
#' This function optimizes the likelihood function using the specified settings.
#' In particular, it implements the the (Nelder-Mead + COBYLA)-hybrid algorithm
#' as described in D'Haen et al. (2025).
#'
#' @param lag_degs Laguerre polynomial maximal degrees.
#' @param cop_name Name of considered copula.
#' @param Y Vector of observed times.
#' @param X Design matrix.
#' @param Delta Vector of censoring indicators.
#' @param init_val_list List of initial values.
#' @param opt_algo Algorithm to be used for optimization. Could be one of
#' \code{"NM"} or \code{"NMCob"}, corresponding to Nelder-Mead and the
#' (Nelder-Mead + COBYLA)-hybrid, respectively.
#' @param hp List of hyperparameters.
#' @param indep_assumption Boolean indicator whether independence can be
#' assumed. Default is \code{indep_assumption = FALSE}.
#'
#' @returns List of program value and argument minimizers (i.e. ML estimates).
#'
#' @noRd
#'
perform_mll_optimisation <- function(lag_degs, cop_name, Y, X, Delta,
                                     init_val_list,
                                     opt_algo, hp,
                                     indep_assumption = FALSE) {

  # Set some parameters
  it_num_unconstr <- hp[["it_NMCob_unconstr"]]
  it_num_constr <- hp[["it_NMCob_constr"]]

  # Set initial values for max/argmax log-likelihood.
  current_mllh <- -Inf
  current_paras <- init_val_list[[1]]

  # Starting from each initial value, optimize the log-likelihood using the
  # Nelder-Mead algorithm.
  if (opt_algo == "NM") {

    # If independence is assumed between T and C...
    if (indep_assumption) {

      # In the case where independence between T and C can be assumed, there is
      # no dependence parameter (first entry of argmax vector) to be estimated.
      current_paras <- current_paras[-1]

      # Optimize the likelihood for each starting value.
      for (init_vals in init_val_list) {
        res <- optim(init_vals[-1], log_lik_fun_indep, method = "Nelder-Mead",
                     control = list(fnscale = -1, maxit = it_num_unconstr),
                     lag_degs = lag_degs, cop_name = cop_name,
                     Y = Y, X = X, Delta = Delta, hp = hp)
        new_paras <- res[[1]]
        new_mllh <- res[[2]]

        # If necessary, update max and argmax.
        if (new_mllh > current_mllh){
          current_mllh <- new_mllh
          current_paras <- new_paras
        }
      }

      # If independence is not assumed between T and C...
    } else {

      # Optimize the likelihood for each starting value.
      for (init_vals in init_val_list) {
        res <- optim(init_vals, log_lik_fun, method = "Nelder-Mead",
                     control = list(fnscale = -1, maxit = it_num_unconstr),
                     lag_degs = lag_degs, cop_name = cop_name,
                     Y = Y, X = X, Delta = Delta, hp = hp)
        new_paras <- res[[1]]
        new_mllh <- res[[2]]

        # If necessary, update max and argmax.
        if (new_mllh > current_mllh){
          current_mllh <- new_mllh
          current_paras <- new_paras
        }
      }
    }
  }

  # Starting from each initial value, optimize the log-likelihood using the
  # Nelder-Mead -- COBYLA hybrid algorithm.
  if (opt_algo == "NMCob"){

    # Set optimization hyperparameters.
    opts <- list("algorithm"= "NLOPT_LN_COBYLA", "maxeval"= it_num_constr)

    if (indep_assumption) {

      # In the case where independence between T and C can be assumed, there is
      # no dependence parameter (first entry of argmax vector) to be estimated.
      current_paras <- current_paras[-1]

      # Optimize the likelihood for each starting value.
      for (init_vals in init_val_list) {

        # Nelder-Mead first stage
        intermed_para <- optim(init_vals[-1], log_lik_fun_indep,
                               method = "Nelder-Mead",
                               control = list(fnscale = -1,
                                              maxit = it_num_unconstr),
                               lag_degs = lag_degs, cop_name = cop_name,
                               Y = Y, X = X, Delta = Delta, hp = hp)$par

        # COBYLA second stage
        res <- nloptr(x0 = intermed_para,
                      eval_f = neg_log_lik_fun_indep, # max. rather than min.!
                      eval_g_ineq = ctuity_constr_fun_indep,
                      opts = opts,
                      lag_degs = lag_degs, cop_name = cop_name,
                      Y = Y, X = X, Delta = Delta, hp = hp)
        new_paras <- res$solution
        new_mllh <- - res$objective # Compare with llh rather than neg_llh

        # Update max/argmax if necessary
        if (new_mllh > current_mllh){
          current_mllh <- new_mllh
          current_paras <- new_paras
        }
      }
    } else {
      for (init_vals in init_val_list){
        intermed_para <- optim(init_vals, log_lik_fun,
                               method = "Nelder-Mead",
                               control = list(fnscale = -1,
                                              maxit = it_num_unconstr),
                               lag_degs = lag_degs, cop_name = cop_name,
                               Y = Y, X = X, Delta = Delta, hp = hp)$par

        res <- nloptr(x0 = intermed_para,
                      eval_f = neg_log_lik_fun, # max. rather than min.!
                      eval_g_ineq = ctuity_constr_fun,
                      opts = opts,
                      lag_degs = lag_degs, cop_name = cop_name,
                      Y = Y, X = X, Delta = Delta, hp = hp)
        new_paras <- res$solution
        new_mllh <- - res$objective # Compare with llh rather than neg_llh

        if (new_mllh > current_mllh){
          current_mllh <- new_mllh
          current_paras <- new_paras
        }
      }
    }
  }

  return(list(current_mllh, current_paras))
}


#' @title Log-likelihood function with subvector of fixed parameters.
#'
#' @description
#' This function is a wrapper around \code{log_lik_fun} implementing the case
#' where a part of the parameter vector is fixed and hence excluded from
#' optimization.
#'
#' @param para_variable Parameters to be estimated.
#' @param para_fixed Fixed parameters.
#' @param fixed_indices Indices of fixed parameters.
#' @param lag_degs Maximal degrees of Laguerre polynomials.
#' @param cop_name Name of considered copula.
#' @param Y Vector of observed times.
#' @param X Design matrix.
#' @param Delta Vector of censoring indicators.
#' @param hp List of hyperparameters.
#'
#' @returns Log-likelihood function evaluation.
#'
#' @seealso log_lik_fun
#'
#' @noRd
#'
log_lik_partial <- function(para_variable, para_fixed, fixed_indices,
                            lag_degs, cop_name, Y, X, Delta, hp) {
  all_paras <- rep(0, length(para_variable) + length(para_fixed))
  all_paras[fixed_indices] <- para_fixed
  all_paras[-fixed_indices] <- para_variable
  return(log_lik_fun(all_paras, lag_degs, cop_name, Y, X, Delta, hp))
}

#' @title Negative of partial log-likelihood function with subvector of fixed
#' parameters.
#'
#' @description
#' This function is a wrapper around \code{log_lik_partial} implementing the case
#' where a part of the parameter vector is fixed and hence excluded from
#' optimization.
#'
#' @inheritParams log_lik_partial
#' @returns Negative log-likelihood function evaluation.
#' @seealso log_lik_partial
#'
#' @noRd
#'
neg_log_lik_partial <- function(para_variable, para_fixed, fixed_indices,
                                lag_degs, cop_name, Y, X, Delta, hp) {
  -log_lik_partial(para_variable, para_fixed, fixed_indices, lag_degs, cop_name,
                   Y, X, Delta, hp)
}

#' @title Continuity constraint of EAL distribution (some fixed parameters).
#'
#' @description
#' This function implements the equality constraint (6.1) of D'Haen et al.
#' (2025) in the case where some parameters are fixed rather than estimated.
#'
#' @param para_variable Parameters to be estimated.
#' @param para_fixed Fixed parameters.
#' @param fixed_indices Indices of fixed parameters.
#' @param lag_degs Maximal degrees of Laguerre polynomials.
#' @param cop_name Name of considered copula.
#' @param Y Vector of observed times.
#' @param X Design matrix.
#' @param Delta Vector of censoring indicators.
#' @param hp List of hyperparameters.
#'
#' @returns Vector of continuity constraint violations.
#'
#' @seealso ctuity_constr_fun
#'
#' @noRd
#'
ctuity_constr_partial <- function(para_variable, para_fixed, fixed_indices,
                                  lag_degs, cop_name, Y, X, Delta, hp) {

  # Extract hyperparameters
  phi_til_index <- hp[["phi_til_index"]]
  phi_index <- hp[["phi_index"]]

  if (cop_name == "indep"){
    ktau_length <- 0
  } else {
    ktau_length <- 1
  }

  all_paras <- rep(0, length(para_variable) + length(para_fixed))
  all_paras[fixed_indices] <- para_fixed
  all_paras[-fixed_indices] <- para_variable

  phis <- all_paras[c(phi_til_index(lag_degs), phi_index(lag_degs))
                    + ktau_length]

  if (lag_degs[1] == 0){
    phi_til_ext <- c(1)
    phi_ext <- c(1, phis)
  } else {
    if (lag_degs[2] == 0){
      phi_til_ext <- c(1, phis)
      phi_ext <- c(1)
    } else {
      phi_til_ext <- c(1, phis[1:lag_degs[1]])
      phi_ext <- c(1, phis[(lag_degs[1]+1):(lag_degs[1] + lag_degs[2])])
    }
  }

  # combination of Lag. polynomials evaluated in 0 on neg./pos. axis:
  left <- sum(phi_til_ext)^2 / sqnorm(phi_til_ext)
  right <- sum(phi_ext)^2 / sqnorm(phi_ext)

  jump <- left - right

  return(c(jump, - jump))
}

#' @title Optimize the log-likelihood function.
#'
#' @description
#' This function optimizes the likelihood function using the specified settings.
#' In particular, it implements the the (Nelder-Mead + COBYLA)-hybrid algorithm
#' as described in D'Haen et al. (2025).
#'
#' @param lag_degs Laguerre polynomial maximal degrees.
#' @param cop_name Name of considered copula.
#' @param Y Vector of observed times.
#' @param X Design matrix.
#' @param Delta Vector of censoring indicators.
#' @param para_fixed Fixed parameters.
#' @param fixed_indices Indices of fixed parameters.
#' @param init_values_variable Initial values for variable parameters.
#' @param opt_algo Algorithm to be used for optimization. Could be one of
#' \code{"NM"} or \code{"NMCob"}, corresponding to Nelder-Mead and the
#' (Nelder-Mead + COBYLA)-hybrid, respectively.
#' @param hp List of hyperparameters.
#' @param indep_assumption Boolean indicator whether independence can be
#' assumed. Default is \code{indep_assumption = FALSE}.
#'
#' @returns List of program value and argument minimizers (i.e. ML estimates).
#'
#' @noRd
#'
perform_mll_partial_optimisation <- function(lag_degs, cop_name, Y, X, Delta,
                                             para_fixed, fixed_indices,
                                             init_values_variable,
                                             opt_algo, hp,
                                             indep_assumption = FALSE) {

  it_num_constr <- hp[["it_NMCob_constr"]]
  it_num_unconstr <- hp[["it_NMCob_unconstr"]]

  current_mllh <- -100000
  current_para_variable <- init_values_variable[[1]]

  if (opt_algo == "NM") {
    if (indep_assumption) {
      current_para_variable <- current_para_variable[-1]
      for (init_vals in init_values_variable){
        res <- optim(init_vals[-1], log_lik_indep_partial, method = "Nelder-Mead",
                     control = list(fnscale = -1, maxit = it_num_unconstr),
                     para_fixed = para_fixed, fixed_indices = fixed_indices,
                     lag_degs = lag_degs, cop_name = cop_name,
                     Y = Y, X = X, Delta = Delta, hp = hp)
        new_para_variable <- res[[1]]
        new_mllh <- res[[2]]

        if (new_mllh > current_mllh){
          current_mllh <- new_mllh
          current_para_variable <- new_para_variable
        }
      }
    } else {
      for (init_vals in init_values_variable){
        res <- optim(init_vals, log_lik_partial, method = "Nelder-Mead",
                     control = list(fnscale = -1, maxit = it_num_unconstr),
                     para_fixed = para_fixed, fixed_indices = fixed_indices,
                     lag_degs = lag_degs, cop_name = cop_name,
                     Y = Y, X = X, Delta = Delta, hp)
        new_para_variable <- res[[1]]
        new_mllh <- res[[2]]

        if (new_mllh > current_mllh){
          current_mllh <- new_mllh
          current_para_variable <- new_para_variable
        }
      }
    }
  }

  if (opt_algo == "NMCob"){

    opts <- list("algorithm"= "NLOPT_LN_COBYLA", "maxeval"= it_num_constr)

    if (indep_assumption) {
      current_para_variable <- current_para_variable[-1]
      for (init_vals in init_values_variable){
        intermed_para_variable <- optim(init_vals[-1], log_lik_indep_partial,
                                        method = "Nelder-Mead",
                                        control = list(fnscale = -1,
                                                       maxit = it_num_unconstr),
                                        para_fixed = para_fixed,
                                        fixed_indices = fixed_indices,
                                        lag_degs = lag_degs,
                                        cop_name = cop_name,
                                        Y = Y, X = X, Delta = Delta, hp = hp)$par

        res <- nloptr(x0 = intermed_para_variable,
                      eval_f = neg_log_lik_indep_partial, # max. rather than min.!
                      eval_g_ineq = ctuity_constr_partial,
                      opts = opts,
                      para_fixed = para_fixed, fixed_indices = fixed_indices,
                      lag_degs = lag_degs, cop_name = cop_name,
                      Y = Y, X = X, Delta = Delta, hp = hp)
        new_para_variable <- res$solution
        new_mllh <- - res$objective # Compare with llh rather than neg_llh

        if (new_mllh > current_mllh){
          current_mllh <- new_mllh
          current_para_variable <- new_para_variable
        }
      }
    } else {
      for (init_vals in init_values_variable){
        intermed_para_variable <- optim(init_vals, log_lik_partial,
                                        method = "Nelder-Mead",
                                        control = list(fnscale = -1,
                                                       maxit = it_num_unconstr),
                                        para_fixed = para_fixed,
                                        fixed_indices = fixed_indices,
                                        lag_degs = lag_degs,
                                        cop_name = cop_name,
                                        Y = Y, X = X, Delta = Delta, hp = hp)$par

        res <- nloptr(x0 = intermed_para_variable,
                      eval_f = neg_log_lik_partial, # max. rather than min.!
                      eval_g_ineq = ctuity_constr_partial,
                      opts = opts,
                      para_fixed = para_fixed, fixed_indices = fixed_indices,
                      lag_degs = lag_degs, cop_name = cop_name,
                      Y = Y, X = X, Delta = Delta, hp = hp)
        new_para_variable <- res$solution
        new_mllh <- - res$objective # Compare with llh rather than neg_llh

        if (new_mllh > current_mllh){
          current_mllh <- new_mllh
          current_para_variable <- new_para_variable
        }
      }
    }
  }

  # Final step: combine results in one vector
  para_all <- rep(0, length(current_para_variable) + length(para_fixed))
  para_all[fixed_indices] <- para_fixed
  para_all[-fixed_indices] <- current_para_variable

  return(list(current_mllh, para_all))
}



#' @title Log-likelihood under independence and with a fixed subvector.
#'
#' @description
#' This function defines the likelihood under the independence assumption and
#' with some elements of the covariate vector being fixed. Specifically, note
#' that \code{para_variable} nor \code{para_fixed} has an entry for \eqn{\eta},
#' the dependence coefficient.
#'
#' @param para_variable Parameters to be estimated.
#' @param para_fixed Fixed parameters.
#' @param fixed_indices Indices of fixed parameters.
#' @param lag_degs Maximal degrees of Laguerre polynomials.
#' @param cop_name Name of considered copula.
#' @param Y Vector of observed times.
#' @param X Design matrix.
#' @param Delta Vector of censoring indicators.
#' @param hp List of hyperparameters.
#'
#' @noRd
#'
log_lik_indep_partial <- function(para_variable, para_fixed, fixed_indices,
                                  lag_degs, cop_name, Y, X, Delta, hp) {
  all_paras <- rep(0, length(para_variable) + length(para_fixed))
  all_paras[fixed_indices] <- para_fixed
  all_paras[-fixed_indices] <- para_variable
  return(log_lik_fun_indep(all_paras, lag_degs, cop_name, Y, X, Delta, hp))
}

#' @title Negative log-likelihood under independence and with a fixed subvector.
#'
#' @description
#' This function defines the negative of the log-likelihood under the
#' independence assumption and with some elements of the covariate vector being
#' fixed. This function is simply a wrapper around \code{log_lik_indep_partial}.
#'
#' @inheritParams log_lik_indep_partial
#' @seealso [log_lik_indep_partial()]
#'
#' @noRd
#'
neg_log_lik_indep_partial <- function(para_variable, para_fixed, fixed_indices,
                                      lag_degs, cop_name, Y, X, Delta, hp) {
  -log_lik_indep_partial(para_variable, para_fixed, fixed_indices, lag_degs,
                         cop_name, Y, X, Delta, hp)
}


#' @title Log-likelihood function under independence.
#'
#' @description
#' This function defines the log-likelihood for all parameters under the
#' independence assumption. Specifically, note that \code{para_variable} nor
#' \code{para_fixed} has an entry for \eqn{\eta}, the dependence coefficient.
#'
#' @param all_paras_indep Vector of parameters (without an entry for \eqn{\eta}
#' in the first position).
#' @param lag_degs Maximal degrees of Laguerre polynomials.
#' @param cop_name Name of considered copula.
#' @param Y Vector of observed times.
#' @param X Design matrix.
#' @param Delta Vector of censoring indicators.
#' @param hp List of hyperparameters.
#'
#' @noRd
#'
log_lik_fun_indep <- function(all_paras_indep, lag_degs, cop_name, Y, X, Delta,
                              hp) {

  #### 0) Extract hyperparameters ####

  beta_full_index <- hp[["beta_full_index"]]
  gamma_index <- hp[["gamma_index"]]
  lambda_index <- hp[["lambda_index"]]
  phi_til_index <- hp[["phi_til_index"]]
  phi_index <- hp[["phi_index"]]
  lambda <- hp[["lambda"]]
  homoscedastic <- hp[["homoscedastic"]]
  heteroscedastic <- !homoscedastic

  #### 1) Retrieving current parameters ####

  paras_T_indices <- c(beta_full_index, gamma_index, lambda_index,
                       phi_til_index(lag_degs), phi_index(lag_degs))

  # no eta present, so no shift + 1
  transf_paras_T <- all_paras_indep[paras_T_indices]
  transf_paras_C <- all_paras_indep[-paras_T_indices]

  beta_full <- transf_paras_T[beta_full_index]
  if (heteroscedastic) {
    transf_lam <- transf_paras_T[lambda_index]
    lambda <- lambda_inverse_transform(transf_lam)
  }
  gam_par <- transf_paras_T[gamma_index]
  phi_til <- transf_paras_T[phi_til_index(lag_degs)]
  phi <- transf_paras_T[phi_index(lag_degs)]

  #### 2) Construct argument matrices for censored & uncensored cases ####

  Y_un <- Y[Delta == 1]
  Y_cens <- Y[Delta == 0]

  X_un <- X[Delta == 1, , drop = FALSE]
  X_cens <- X[Delta == 0, , drop = FALSE]

  #### 3) Density & distribution evaluation ####

  # density
  dens_T_un <- dEAL(Y_un, X_un, beta_full, gam_par, lambda, phi_til, phi, hp)
  dens_C_cens <- dCens(Y_cens, X_cens, transf_paras_C)

  # distribution
  distr_T_cens <- pEAL(Y_cens, X_cens, beta_full, gam_par, lambda, phi_til, phi)
  distr_C_un <- pCens(Y_un, X_un, transf_paras_C)

  #### 4) Compute and return log-likelihood ####

  lh_uncens <- dens_T_un * (1 - distr_C_un)
  lh_cens <- dens_C_cens * (1 - distr_T_cens)

  llh <- sum(log(lh_uncens)) + sum(log(lh_cens))

  res <- min(max(llh, -100000), 100000) # avoid returning +- Inf
  if (is.nan(res)){ # avoid NaN (due to Inf - Inf by num. errors) as well
    res <- -100000
  }
  return(res)
}

#' @title Negative log-likelihood function under independence.
#'
#' @description
#' This function defines the negative log-likelihood for all parameters under
#' the independence assumption. Specifically, note that \code{para_variable} nor
#' \code{para_fixed} has an entry for \eqn{\eta}, the dependence coefficient.
#' This function is simply a wrapper around \code{log_lik_fun_indep}.
#'
#' @inheritParams log_lik_fun_indep
#' @seealso [neg_log_lik_fun_indep()]
#'
#' @noRd
#'
neg_log_lik_fun_indep <- function(all_paras_indep, lag_degs, cop_name, Y, X,
                                  Delta, hp) {
  -log_lik_fun_indep(all_paras_indep, lag_degs, cop_name, Y, X, Delta, hp)
}

#' @title Continuity constraint of EAL distribution under independence.
#'
#' @description
#' Functionality is similar to \code{ctuity_constr_fun}.
#'
#' @param all_paras_indep Parameters to be estimated.
#' @param lag_degs Maximal degrees of Laguerre polynomials.
#' @param cop_name Name of considered copula.
#' @param Y Vector of observed times.
#' @param X Design matrix.
#' @param Delta Vector of censoring indicators.
#' @param hp List of hyperparameters.
#'
#' @returns Vector of continuity constraint violations.
#'
#' @seealso ctuity_constr_fun
#'
#' @noRd
#'
ctuity_constr_fun_indep <- function(all_paras_indep, lag_degs, cop_name, Y, X,
                                    Delta, hp) {

  # Extract hyperparameters
  phi_til_index <- hp[["phi_til_index"]]
  phi_index <- hp[["phi_index"]]

  phi_til <- all_paras_indep[phi_til_index(lag_degs)]
  phi <- all_paras_indep[phi_index(lag_degs)]

  phi_til_ext <- c(1, phi_til)
  phi_ext <- c(1, phi)

  # combination of Lag. polynomials evaluated in 0 on neg./pos. axis:
  left <- sum(phi_til_ext)^2 / sqnorm(phi_til_ext)
  right <- sum(phi_ext)^2 / sqnorm(phi_ext)

  jump <- left - right

  return(c(jump, - jump))
}


#' @title Degree selection based on AIC.
#'
#' @description
#' This function uses AIC to select the degree of the Laguerre polynomials in
#' the enriched Laplace distribution.
#'
#' @param Y Vector of observed times.
#' @param X matrix of observed covariate values.
#' @param T1 Subvector of observed times corresponding to \eqn{Delta = 1}, i.e.
#' vector of observed event times.
#' @param X1 Subvector of covariates corresponding to \eqn{Delta = 1}.
#' @param C0 Subvector of observed times corresponding to \eqn{Delta = 0}, i.e.
#' vector of observed censoring times.
#' @param X0 Subvector of covariates corresponding to \eqn{Delta = 1}.
#' @param Delta Vector of censoring indicators.
#' @param cop_name Name of the copula to use.
#' @param hp List of hyperparameters.
#'
#' @noRd
#'
deg_sel_AIC_grid_only_T <- function(cop_name, Y, X, Delta,
                                    T1, X1, C0, X0, hp,
                                    indep_assumption = FALSE) {

  #### 0) Extract hyperparameters ####

  sample_size <- hp[["sample_size"]]
  max_deg <- hp[["lag_deg_upper"]]
  init_value_number_intermed <- hp[["init_value_number_intermed"]]
  init_value_number_basis <- hp[["init_value_number_basis"]]
  lower_bound_on_phis <- hp[["lower_bound_on_phis"]]
  upper_bound_on_phis <- hp[["upper_bound_on_phis"]]
  it_NM <- hp[["it_NM"]]
  beta_full_index <- hp[["beta_full_index"]]
  gamma_index <- hp[["gamma_index"]]
  lambda_index <- hp[["lambda_index"]]
  ktau_length <- 1 - as.numeric(indep_assumption)

  #### 1) first perform optimization under zero degrees (using Nelder-Mead) ####

  init_val_list <- determine_init_basis(T1, X1, C0, X0, Delta, cop_name, hp)
  mll_basis <- perform_mll_optimisation(c(0, 0), cop_name,
                                        Y, X, Delta,
                                        init_val_list = init_val_list,
                                        opt_algo = "NM", hp,
                                        indep_assumption = indep_assumption)
  para_basis <- mll_basis[[2]]

  # Benchmark: AIC/BIC score for mll optimisation under degrees (0,0)
  basis_para_num <- length(para_basis)
  BIC_factor <- log(sample_size)
  basis_AIC <- 2*basis_para_num - 2* mll_basis[[1]]
  basis_BIC <- BIC_factor*basis_para_num - 2* mll_basis[[1]]

  T_indices <- ktau_length + c(beta_full_index, gamma_index, lambda_index)
  para_T_basis <- para_basis[T_indices]
  para_fixed_basis <- para_basis[-T_indices]

  # AIC/BIC degree selection
  all_AIC_scores <- matrix(rep(0, (max_deg + 1)^2),
                           nrow = (max_deg + 1), byrow = TRUE)
  all_AIC_scores[1,1] <- basis_AIC
  best_AIC_degrees <- c(0,0)
  best_AIC_score <- basis_AIC
  best_AIC_transf_paras <- para_basis

  all_BIC_scores <- matrix(rep(0, (max_deg + 1)^2),
                           nrow = (max_deg + 1), byrow = TRUE)
  all_BIC_scores[1,1] <- basis_BIC
  best_BIC_degrees <- c(0,0)
  best_BIC_score <- basis_BIC
  best_BIC_transf_paras <- para_basis

  #### 2) Then perform optimisation under range of degrees (grid approach) ####

  for (L_degree in (seq_len(max_deg + 1) - 1)) {
    for (R_degree in (seq_len(max_deg + 1) - 1)) {

      if (L_degree + R_degree > 0) {

        # Determine initial values
        total_degree <- L_degree + R_degree

        # fix everything except the indices for T (new phi's also not fixed)
        fixed_indices <- seq_len(
          length(para_basis) + total_degree
        )[-c(T_indices, max(T_indices) + seq_len(total_degree))]

        init_val_variable <- vector(mode = "list",
                                    length = init_value_number_intermed)
        for (j in 1:init_value_number_intermed){
          random_phis <- runif(total_degree,
                               lower_bound_on_phis, upper_bound_on_phis)

          if (ktau_length == 1) {
            init_val_variable[[j]] <- c(
              slight_perturb(para_T_basis),
              random_phis)
          } else {
            init_val_variable[[j]] <- c(0, slight_perturb(para_T_basis),
                                        random_phis)
          }
        }

        new_mll <- perform_mll_partial_optimisation(
          lag_degs = c(L_degree, R_degree),
          cop_name = cop_name,
          Y = Y, X = X, Delta = Delta,
          para_fixed = para_fixed_basis,
          fixed_indices = fixed_indices,
          init_values_variable = init_val_variable,
          opt_algo = "NMCob",
          hp = hp,
          indep_assumption = indep_assumption)

        new_AIC <- 2*(basis_para_num + total_degree) - 2* new_mll[[1]]
        new_BIC <- BIC_factor*(basis_para_num + total_degree) - 2* new_mll[[1]]

        # Storing result in grid
        all_AIC_scores[L_degree + 1, R_degree + 1] <- new_AIC
        all_BIC_scores[L_degree + 1, R_degree + 1] <- new_BIC

        # Update minimum
        if (new_AIC < best_AIC_score){
          best_AIC_degrees <- c(L_degree, R_degree)
          best_AIC_score <- new_AIC
          best_AIC_transf_paras <- new_mll[[2]]
        }
        if (new_BIC < best_BIC_score){
          best_BIC_degrees <- c(L_degree, R_degree)
          best_BIC_score <- new_BIC
          best_BIC_transf_paras <- new_mll[[2]]
        }
      }
    }
  }

  #### 3) Return the results ####

  return(list(best_AIC_degrees, best_AIC_transf_paras))
}


#' @title Validate input to the function \code{QRdepCens}.
#'
#' @description
#' This function checks a set of preconditions that are required for the
#' function \code{QRdepCens} to run successfully.
#'
#' @param data Data frame.
#' @param hp List of hyperparameters. See the documentation of \code{QRdepCens}
#' for an exhaustive summary of options that can be specified.
#'
#' @noRd
#'
QRdepCens_validateInput <- function(data, hp) {

  #### Checks related to the data ####

  # Initialize vector of column indices which have not been checked.
  cols.to.check <- 1:ncol(data)

  # The data should contain a column pertaining to the observed times, named 'Y'
  if ("Y" %notin% colnames(data)) {
    stop("The column pertaining to the observed times should be named 'Y'.")
  }
  cols.to.check <- setdiff(cols.to.check, which(colnames(data) == "Y"))

  # The data should contain a column pertaining to the censoring indicator,
  # named Delta. Values in this column can only be 0 or 1.
  if ("Delta" %notin% colnames(data)) {
    stop("The column pertaining to the censoring indicator should be named 'Delta'.")
  } else if (any(data$Delta %notin% c(0, 1))) {
    stop("The censoring indicator should be either 0 or 1.")
  }
  cols.to.check <- setdiff(cols.to.check, which(colnames(data) == "Delta"))

  # The data should contain at least one covariate. Columns pertaining to
  # covariates should be named X1, X2, etc.
  if ("X1" %notin% colnames(data)) {
    stop("The data should contain at least one covariate, which should be named 'X1'.")
  }
  cols.to.check <- setdiff(cols.to.check, which(grepl("X[1-9][[:digit:]]*$", colnames(data))))

  # The data should not already contain an intercept column
  if (any(apply(extract.covariates(data), 2, function(col) {all(col == 1)}))) {
    stop("The data should not already contain an intercept column.")
  }

  # The data should not contain any other columns
  if (length(cols.to.check) != 0) {
    stop("Unrecognized columns detected. See documentation of 'QRdepCens'.")
  }

  #### Checks related to the hyperparameters ####

  # The copula to be used should be specified, and be one of the available
  # options.
  if (is.null(hp[["test_cop_name"]])) {
    stop("The copula to be used should be specified.")
  } else if (hp[["test_cop_name"]] %notin% c("frank", "frankPos", "gumbel", "clayton", "indep")) {
    stop("Specified copula not implemented. See documentation of QRdepCens.")
  }

  # The maximal degree of the Laguerre polynomials should specified, and be a
  # positive integer.
  if (is.null(hp[["lag_deg_upper"]])) {
    stop("The maximum degree of Laguerre polynomials should be specified.")
  } else if (hp[["lag_deg_upper"]] <= 0) {
    stop("Maximum degree of Laguerre polynomial should be positive.")
  } else if (!is_integer(hp[["lag_deg_upper"]])) {
    stop("Maximum degree of Laguerre polynomial should be integer.")
  }

  # The maximum number of contraint optimization iterations should be specified
  # and positive.
  if (is.null(hp[["it_NMCob_constr"]])) {
    stop("The maximum number of contraint optimization iterations should be specified.")
  } else if (hp[["it_NMCob_constr"]] <= 0) {
    stop("The maximum number of contraint optimization iterations should be positive.")
  }

  # The maximum number of uncontraint optimization iterations should be
  # specified and positive.
  if (is.null(hp[["it_NMCob_unconstr"]])) {
    stop("The maximum number of uncontraint optimization iterations should be specified.")
  } else if (hp[["it_NMCob_unconstr"]] <= 0) {
    stop("The maximum number of uncontraint optimization iterations should be positive.")
  }

  # The value indicating if homoscedasticity can be assumed should be boolean.
  if (!is.logical(hp[["homoscedastic"]])) {
    stop("'homoscedastic' should be logical.")
  }

  # Check whether the quantile of interest/quantile index has been correctly
  # specified (depending on whether homoscedasticity is assumed)
  if (hp[["homoscedastic"]]) {
    if (!is.null(hp[["lambda_index"]])) {
      stop("'lambda_index' should not be specified for a homoscedastic model.")
    }
    if (is.null(hp[["lambda"]])) {
      stop("quantile of interest 'lambda' should be specified in homoscedastic model.")
    } else if (!is.numeric(hp[["lambda"]])) {
      stop("'lambda' should be numeric.")
    } else if (hp[["lambda"]] < 0 | hp[["lambda"]] > 1) {
      stop("'lambda' should be a value between 0 and 1.")
    }
  } else {
    if (is.null(hp[["lambda_index"]])) {
      stop("'lambda_index' should be specified for heteroscedastic model.")
    }
    if (!is.null(hp[["lambda"]])) {
      stop("'lambda' should not be specified for heteroscedastic model.")
    }
  }

  # Number of variance bootstrap iterations should be a positive integer.
  if (!is.null(hp[["variance.bootstrap.iterations"]])) {
    if (!is_integer(hp[["variance.bootstrap.iterations"]])) {
      stop("The specified number of variance bootstrap iterations should be integer.")
    } else if (hp[["variance.bootstrap.iterations"]] <= 1) {
      stop("The specified number of variance bootstrap iterations should be strictly larger than one.")
    }
  }

  # ToDo: init_value_number_basis, variance.b

}


#' @title Obtain default settings for hyperparameters.
#'
#' @description
#' This function returns the default settings for the hyperparameters used in
#' running the copula-based quantile regression method.
#'
#' @param data Data frame.
#' @param homoscedastic Boolean value indicating whether homoscedasticity can
#' be assumed for the model. Default is \code{homoscedastic = FALSE}.
#'
#' @returns A list containing the elements:
#' \describe{
#'  \item{homoscedastic}{Boolean value indicating whether homoscedasticity can
#'  be assumed for the model.}
#'  \item{lambda}{Quantile of interest.}
#'  \item{test_cop_name}{Name of copula to be used.}
#'  \item{lag_deg_upper}{Upper bound on Laguerre polynomial degree to be
#'  considered.}
#'  \item{it_NMCob_constr, it_NMCob_unconstr, it_NM:}{Maximum number of
#'  optimization cycles to run during runtime.}
#'  \item{init_value_number_basis, init_value_number_intermed,
#'  init_value_number_final:}{Number of initial values to be considered in
#'  optimization routines.}
#'  \item{lower_bound_on_phis, upper_bound_on_phis}{Lower and upper bound on
#'  the initial values for the Laguerre polynomial coefficients.}
#'  \item{variance.bootstrap.iterations:}{Number of bootstrap iterations during
#'  variance estimation.}
#'  \item{sample_size:}{The sample size.}
#'  \item{X_dimension:}{The number of covariates in the model.}
#'  \item{beta_full_index:}{Vector of indices in the full parameter vector
#'  containing the elements of beta.}
#'  \item{gamma_index:}{Vector of indices in the full parameter vector
#'  containing the element(s) (of) gamma.}
#'  \item{lambda_index:}{Index in the full parameter vector for quantile of
#'  interest.}
#'  \item{phi_til_index, phi_index:}{Index in the full parameter vector for
#'  Laguerre polynomial coefficients.}
#' }
#'
#' @noRd
#'
QRdepCens_getDefaultHyperparameters <- function(data, homoscedastic = FALSE) {

  # Set default hyperparameters
  test_cop_name <- "frank"
  lag_deg_upper <- 2
  it_NMCob_constr <- 100
  it_NMCob_unconstr <- 100
  it_NM <- 100
  init_value_number_basis <- 5
  init_value_number_intermed <- 5
  init_value_number_final <- 5
  lower_bound_on_phis <- -2
  upper_bound_on_phis <- 2
  variance.bootstrap.iterations <- 50

  # Define the positions of elements of subvectors of the parameter vector that
  # pertain to certain parts of the model for T.
  sample_size <- nrow(data)
  X_dimension <- ncol(extract.covariates(data))
  beta_full_index <- 1:(X_dimension + 1)

  # Depending on whether or not homoscedasticity is assumed, set the indices
  # of the parameters in the full parameter vector.
  if (homoscedastic) {
    gamma_index <- X_dimension + 2
    lambda <- 0.5
    lambda_index <- NULL
    phi_til_index <- function(lag_degs) {
      if (lag_degs[1] == 0){
        return(c())
      } else {
        return((X_dimension + 3):(X_dimension + 2 + lag_degs[1]))
      }
    }
    phi_index <- function(lag_degs) {
      if (lag_degs[2] == 0){
        return(c())
      } else {
        return((X_dimension + 3 + lag_degs[1]):
                 (X_dimension + 2 + lag_degs[1] + lag_degs[2]))
      }
    }
  } else {
    gamma_index <- (X_dimension + 2) : (2*X_dimension + 2)
    lambda <- NULL
    lambda_index <- c(2*X_dimension + 3)
    phi_til_index <- function(lag_degs) {
      if (lag_degs[1] == 0){
        return(c())
      } else {
        return((2*X_dimension + 4):(2*X_dimension + 3 + lag_degs[1]))
      }
    }
    phi_index <- function(lag_degs) {
      if (lag_degs[2] == 0){
        return(c())
      } else {
        return((2*X_dimension + 4 + lag_degs[1]):
                 (2*X_dimension + 3 + lag_degs[1] + lag_degs[2]))
      }
    }
  }

  # Return the result
  list(homoscedastic = homoscedastic,
       lambda = lambda,
       test_cop_name = test_cop_name,
       lag_deg_upper = lag_deg_upper,
       it_NMCob_constr = it_NMCob_constr,
       it_NMCob_unconstr = it_NMCob_unconstr,
       it_NM = it_NM,
       init_value_number_basis = init_value_number_basis,
       init_value_number_intermed = init_value_number_intermed,
       init_value_number_final = init_value_number_final,
       lower_bound_on_phis = lower_bound_on_phis,
       upper_bound_on_phis = upper_bound_on_phis,
       variance.bootstrap.iterations = variance.bootstrap.iterations,
       sample_size = sample_size,
       X_dimension = X_dimension,
       beta_full_index = beta_full_index,
       gamma_index = gamma_index,
       lambda_index = lambda_index,
       phi_til_index = phi_til_index,
       phi_index = phi_index)
}

#' @title Print the results of the QRdepCens model.
#'
#' @description
#' This function prints the results of the QRdepCens model to the console.
#'
#' @param results List containing the estimated model parameters, as returned by
#' \code{QRdepCens}.
#'
#' @noRd
#'
QRdepCens_printSummary <- function(results) {

  #### Unpack results ####

  est_beta_full <- results[["est_beta"]]
  est_gam <- results[["est_gamma"]]
  est_lambda <- results[["est_lambda"]]
  cov_mat <- results[["cov_mat"]]

  # Parse variance information
  if (!is.null(cov_mat)) {
    vars <- diag(cov_mat)
    b1 <- 1
    b2 <- length(est_beta_full)
    g1 <- b2 + 1
    g2 <- b2 + length(est_gam)
    l1 <- g2 + 1
    var_beta_full <- vars[b1:b2]
    var_gam <- vars[g1:g2]
    var_lambda <- vars[l1]
  } else {
    var_beta_full <- NA
    var_gam <- NA
    var_lambda <- NA
  }

  #### Arrange all results in a data frame ####

  # Results for beta.
  df.beta <- data.frame(
    "coef" = sprintf("beta%d", 0:(length(est_beta_full) - 1)),
    "est"  = est_beta_full,
    "s.d." = sqrt(var_beta_full),
    "conf. int." = make_ci(est_beta_full, sqrt(var_beta_full))
  )

  # Results for gamma
  df.gam <- data.frame(
    "coef" = sprintf("gamma%d", 0:(length(est_gam) - 1)),
    "est"  = est_gam,
    "s.d." = sqrt(var_gam),
    "conf. int." = make_ci(est_gam, sqrt(var_gam))
  )

  # Results for lambda.
  df.lambda <- data.frame(
    "coef" = "lambda",
    "est"  = est_lambda,
    "s.d." = sqrt(var_lambda),
    "conf. int." = make_ci(est_lambda, sqrt(var_lambda))
  )

  # Concatenate all data frames and convert all columns to string values
  df.res <- rbind(df.beta, df.gam, df.lambda)
  df.res[, 2] <- sprintf("%.3f", df.res[, 2])
  df.res[, 3] <- sprintf("%.3f", df.res[, 3])

  #### Print the results ####

  # Add header
  df.res <- rbind(c("coef", "est", "s.d.", "conf. int."), df.res)

  # For each column, determine the maximal string length
  max.str.lengths <- apply(df.res, 2, function(col) {max(nchar(col))})

  # Set each element of the results data frame to have this column-wise
  # maximal string length.
  for (col.idx in 1:ncol(df.res)) {
    max.length <- max.str.lengths[col.idx]
    sprintf.arg <- paste0("%", max.length, "s")
    df.res[, col.idx] <- sprintf(sprintf.arg, df.res[, col.idx])
  }

  # Define sequence of lines to print
  df.res.str <- apply(df.res, 1, function(row) {sprintf("%s  %s  %s  %s", row[1], row[2], row[3], row[4])})

  # Print the summary
  for (line in df.res.str) {
    message(line)
  }
}

#' @title Estimate the model of D'Haen et al. (2025).
#'
#' @description
#' This function estimates the parameters in the model of D'Haen et al. (2025).
#'
#' @param data Data on which the model should be estimated. Note that the data
#' should be structured in a specific form: The observed times should be put in
#' a column named \code{"Y"}, the censoring indicators in a column named
#' \code{"Delta"}, and the covariates in columns named \code{"X1"}, \code{"X2"},
#' ... (in increasing order). The given data set cannot contain any other
#' columns.
#' @param hp List of hyperparameters to be used, the elements of which will
#' overwrite the default settings. In particular, consider changing:
#' \describe{
#'  \item{test_cop_name:}{Copula to be used. Should be one of
#'  \code{"frank"}, \code{"gumbel"}, \code{"clayton"} or \code{"indep"}.}
#'  \item{homoscedastic:}{Boolean flag indicating whether
#'  homoscedasticity can be assumed.}
#'  \item{variance.bootstrap.iterations:}{Number of bootstrap resamples
#'  to use during variance estimation. Consider increase if more precision is
#'  needed; consider decreasing to reduce computation time.}
#' }
#' Other hyperparameters can be changed though it is not recommended. We refer
#' to the source code for the available options.
#' @param var.estimate Boolean value indicating whether the variance should be
#' estimated (via bootstrap). This can take a considerable amount of time.
#' Default is \code{var.estimate = FALSE}.
#' @param verbose Verbosity flag (boolean) indicating whether the results should
#' be printed to the console. Default is \code{verbose = TRUE}.
#'
#' @note The variance estimation procedure could easily be paralelized.
#' However, this is currently not implemented.
#'
#' @examples
#' \donttest{
#'
#'  # Load the data
#'  data(liver)
#'
#'  # Give standard column names (required!)
#'  colnames(liver) <- c("patient", "Y", "Delta", "X1", "X2", "X3", "X4")
#'  liver <- liver[, c(-1, -6, -7)]
#'
#'  # Run the model
#'  hp <- list(
#'    homoscedastic = FALSE,
#'    test_cop_name = "frank"
#'  )
#'  QRdepCens(liver, hp, var.estimate = FALSE)
#'  # Takes a while if var.estimate = TRUE...
#' }
#'
#' @references D'Haen, M., Van Keilegom, I. and Verhasselt, A. (2025). Quantile
#' regression under dependent censoring with unknown association. Lifetime Data
#' Analysis 31(2):253-299.
#'
#' @export
#'
QRdepCens <- function(data, hp, var.estimate = FALSE, verbose = TRUE) {

  #### Load dependencies ####

  lapply(QRdepCens.getDependencies(), require, character.only = TRUE)

  #### Set hyperparameters ####

  homoscedastic <- if ("homoscedastic" %in% names(hp)) hp$homoscedastic else FALSE
  user.hp <- hp
  hp <- QRdepCens_getDefaultHyperparameters(data, homoscedastic)
  hp[names(user.hp)] <- user.hp

  #### Input validation ####

  QRdepCens_validateInput(data, hp)

  #### Parse data and hyperparameters, and define settings ####

  # Extract different parts of the data
  X <- cbind(1, as.matrix(extract.covariates(data)))
  Y <- data$Y
  Delta <- data$Delta

  # Subvector of observed survival times
  T1 = Y[Delta == 1]
  X1 = X[Delta == 1, , drop = FALSE]

  # Subvector of observed censoring times
  C0 = Y[Delta == 0]
  X0 = X[Delta == 0, , drop = FALSE]

  # Hyperparameters
  cop_name <- hp[["test_cop_name"]]
  lag_deg_upper <- hp[["lag_deg_upper"]]
  it_NMCob_constr <- hp[["it_NMCob_constr"]]
  it_NMCob_unconstr <- hp[["it_NMCob_unconstr"]]
  X_dimension <- hp[["X_dimension"]]

  # Dependence parameter
  if (cop_name == "indep"){
    ktau_length_test <- 0
    indep_assumption <- TRUE
  } else {
    ktau_length_test <- 1
    indep_assumption <- FALSE
  }

  #### Estimate the model ####

  # Determine the optimal Laguerre degrees using AIC
  intermed_res <- deg_sel_AIC_grid_only_T(cop_name = cop_name,
                                          Y = Y, X = X, Delta = Delta,
                                          T1, X1, C0, X0, hp,
                                          indep_assumption = indep_assumption)
  lag_degrees <- intermed_res[[1]]
  current_paras <- intermed_res[[2]]

  # Set maximum likelihood optimization starting values
  final_init_values <- determine_init_final(current_paras, hp, indep_assumption)

  # Estimate the model
  mll_res <- perform_mll_optimisation(lag_degrees, cop_name,
                                      Y, X, Delta,
                                      init_val_list = final_init_values,
                                      opt_algo = "NMCob",
                                      hp = hp,
                                      indep_assumption = indep_assumption)

  # Extract the results
  est_mll <-  mll_res[[1]]
  est_transf_para_vec <- mll_res[[2]]
  est_para_vec <- para_inv_transform(est_transf_para_vec,
                                     cop_name = cop_name,
                                     hp = hp)

  #### Parse the results of model estimation ####

  # Define the positions of elements of subvectors of the parameter vector that
  # pertain to certain parts of the model for T.
  # Define functions that return the indices corresponding to the Laguerre
  # polynomial coefficients in the parameter vector.
  beta_full_index <- hp[["beta_full_index"]]
  gamma_index <- hp[["gamma_index"]]
  lambda_index <- hp[["lambda_index"]]
  phi_til_index <- hp[["phi_til_index"]]
  phi_index <- hp[["phi_index"]]

  # Indices in the full coefficient vector estimate corresponding to the model
  # for T. Note that ktau_length = 0 for when assuming independence between T
  # and C (in which case the dependence parameter was not estimated), and
  # ktau_length = 1 otherwise.
  paras_T_indices <- c(beta_full_index, gamma_index, lambda_index,
                       phi_til_index(lag_degrees),
                       phi_index(lag_degrees)) + ktau_length_test
  est_paras_T <- est_para_vec[paras_T_indices]

  # Extract subvectors of parameters corresponding to different 'parts' of the
  # model for T.
  est_beta_full <- est_paras_T[beta_full_index]
  est_gam_par <- est_paras_T[gamma_index]
  est_lambda <- est_paras_T[lambda_index]
  est_phi_til <- est_paras_T[phi_til_index(lag_degrees)]
  est_phi <- est_paras_T[phi_index(lag_degrees)]

  #### Estimate the variance ####

  # Initialize a NULL value for the covariance matrix.
  cov.mat <- NULL

  if (var.estimate) {

    # Extract number of bootstrap iterations
    B <- hp[["variance.bootstrap.iterations"]]

    # If this number was not specified, throw an error
    if (is.null(B)) {
      stop("Number of variance bootstrap iterations should be specified.")
    }

    # Initialize matrix to store the results
    boot.mat <- NULL

    # Obtain bootstrap estimates of the parameters
    for (b in seq_len(B)) {

      # Construct bootstrapped data set
      data.b <- data[sample(1:nrow(data), nrow(data), replace = TRUE), ]

      # Estimate the parameters
      res.b <- QRdepCens(data.b, hp, var.estimate = FALSE, verbose = FALSE)
      est.b <- c(res.b$est_beta, res.b$est_gamma, res.b$est_lambda)

      # Store the results
      boot.mat <- rbind(boot.mat, est.b)
    }

    # Compute bootstrap covariance matrix
    cov.mat <- var(boot.mat)
  }

  #### Print summary and return the results ####

  # List of results
  rtrn <- list(
    etimate_full = est_para_vec,
    estimate_T_full = est_paras_T,
    est_beta = est_beta_full,
    est_gamma = est_gam_par,
    est_lambda = est_lambda,
    est_phi_til = est_phi_til,
    est_phi = est_phi,
    cov_mat = cov.mat
  )

  # If the homoscedastic option was selected, add the selected quantile of
  # interest to this vector (in place of where its estimate would otherwise go).
  if (hp$homoscedastic) {
    rtrn[["est_lambda"]] <- hp$lambda
    if (var.estimate) {
      rtrn[["cov_mat"]]
    }
  }

  # Print a summary to the console
  if (verbose) {
    QRdepCens_printSummary(rtrn)
  }

  # Return the results (invisibly)
  return(invisible(rtrn))

}

#' @title Estimate quantile based on QRdepCens model.
#'
#' @description
#' This function estimates a specified quantile conditionally on the given
#' covariate values and estimated model parameters. We refer to equation (3.10)
#' in D'haen et al. (2024) for details.
#'
#' @param p Quantile to be estimated.
#' @param x Vector of covariate values at which to estimate the specified
#' quantiles (should not contain an intercept entry).
#' @param res List containing the model parameters, as returned by the function
#' \code{QRdepCens}.
#'
#' @note For the moment, this function is not made visible to the user.
#'
#' @noRd
#'
QRdepCens_estimateQuantile <- function(p, x, res) {

  # Extract the covariates pertaining to T from the model estimates
  beta <- res[["est_beta"]]
  gamma <- res[["est_gamma"]]
  lambda <- res[["est_lambda"]]
  phi_til <- res[["phi_til"]]
  phi <- res[["phi"]]

  # Add intercept to x
  x <- c(1, x)

  # Compute conditional variance of T
  if (length(gamma) == 1) {
    sd_T <- exp(gamma)
  } else {
    sd_T <- exp(x %*% gamma)
  }

  # Compute quantile of error term
  q.error <- qEAL(p, x = c(1), beta = c(0), gam_par = 0, lambda = lambda,
                  phi_til = phi_til, phi = phi)

  # Compute and return quantile estimate
  x %*% beta + sd_T * q.error
}





