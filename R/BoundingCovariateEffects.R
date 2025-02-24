
#### Dependencies ####

# require("survival")
# require("EnvStats")
# require("splines2")
# require("copula")
# require("nloptr")
# require("SPOT")
# require("doParallel")
# require("lubridate")
# require("R6")

#### Chronometer ####

#' @title Chronometer object
#'
#' @description
#' R6 object that mimics a chronometer. It can be started, paused, record legs
#' and stopped.
#'
#' @importFrom lubridate now
#' @importFrom R6 R6Class
#'
#' @export
Chronometer <- R6Class(
  "Chronometer",

  public = list(

    #' @description Display the values stored in this chronometer object.
    show = function() {
      print(paste0("Start: ", private$start.time.abs))
      if (length(private$leg.times.rel) != 0) {
        for (i in 1:length(private$leg.times.rel)) {
          leg.name <- colnames(private$leg.times.rel)[i]
          print(paste0(leg.name, ": ", private$leg.times.rel[i]))
        }
      }
      print(paste0("Stop: ", private$stop.time.rel))
    },

    #' @description Reset the chronometer.
    reset = function() {
      private$start.time.abs = NULL
      private$stop.time.abs = NULL
      private$leg.times.abs = matrix(nrow = 1, ncol = 0)
      private$stop.time.rel = NULL
      private$leg.times.rel = matrix(nrow = 1, ncol = 0)
    },

    #' @description Start the chronometer
    start = function() {
      if (!is.null(private$start.time.abs)) {
        stop("Chronometer already started")
      }
      private$start.time.abs <- lubridate::now()
    },

    #' @description Stop the chronometer
    #' @param leg.name (optional) Name for the stopped leg.
    stop = function(leg.name = NULL) {
      if (is.null(private$start.time.abs)) {
        stop("Chronometer not started yet")
      }
      if (!is.null(private$stop.time.rel)) {
        stop("Chronometer already stopped")
      }

      private$stop.time.abs <- lubridate::now()
      private$stop.time.rel <-
        private$get.time.diff(private$start.time.abs, private$stop.time.abs)

      # If legs were being recorded, also record the final leg time
      if (ncol(private$leg.times.rel) != 0) {
        private$append.leg.time.abs(private$stop.time.abs, leg.name)
        private$update.leg.times.rel()
      }
    },

    #' @description Record a leg time. The chronometer will continue running.
    #' @param leg.name Name for the recorded leg.
    record.leg = function (leg.name = NULL) {
      if (is.null(private$start.time.abs)) {
        stop("Chronometer not started yet")
      }
      if (!is.null(private$stop.time.rel)) {
        stop("Chronometer already stopped")
      }

      # Compute and record leg time
      private$append.leg.time.abs(lubridate::now(), leg.name)
      private$update.leg.times.rel()
    },

    #' @description Like \code{show} method, but more rudimentary.
    get.chronometer.data = function() {
      list(private$start.time.abs,
           private$stop.time.abs,
           private$leg.times.abs,
           private$stop.time.rel,
           private$leg.times.rel
      )
    },

    #' @description Return the total time span between start and stop.
    #' @param force Boolean variable. If \code{TRUE}, avoids error when calling
    #' this function while chronometer has not been stopped yet.
    get.total.time = function(force = FALSE) {
      if (is.null(private$stop.time.rel)) {
        if (force) {
          return(-Inf)
        } else {
          stop("Chronometer not stopped yet")
        }
      }
      private$stop.time.rel
    },

    #' @description Return total time spent per leg category (using leg names).
    #' @param force force Boolean variable. If \code{TRUE}, avoids error when calling
    #' this function while chronometer has not been stopped yet.
    accumulate.legs = function(force = FALSE) {
      if (ncol(private$leg.times.rel) == 0) {
        if (force) {
          return(rep(-Inf, 4))
        } else {
          stop("No leg times recorded yet.")
        }
      }
      categories <- unique(colnames(private$leg.times.rel))
      rtrn <- matrix(nrow = 1, ncol = length(categories))
      colnames(rtrn) <- categories
      for (cat in categories) {
        cat.idx <- which(categories == cat)
        idxs.col.to.sum <- which(colnames(private$leg.times.rel) == cat)
        rtrn[1, cat.idx] <- sum(private$leg.times.rel[1, idxs.col.to.sum])
      }
      rtrn
    }
  ),

  private = list(

    # Variables
    start.time.abs = NULL,
    stop.time.abs = NULL,
    leg.times.abs = matrix(nrow = 1, ncol = 0),
    stop.time.rel = NULL,
    leg.times.rel = matrix(nrow = 1, ncol = 0),

    # Function to get time differences in seconds
    get.time.diff = function(t1, t2) {
      as.numeric(abs(difftime(t1, t2, units = "secs")))
    },

    # Function to append a leg time to the matrix of leg times.
    append.leg.time.abs = function(time, leg.name) {
      old.col.names <- colnames(private$leg.times.abs)
      if (is.null(leg.name)) {
        leg.name <- sprintf("Leg %s", ncol(private$leg.times.abs) + 1)
      }
      private$leg.times.abs <- cbind(private$leg.times.abs, time)
      colnames(private$leg.times.abs) <- c(old.col.names, leg.name)
    },

    # Function to update the relative leg times based on the absolute leg times
    update.leg.times.rel = function () {
      time.vct <- c(private$start.time.abs, c(private$leg.times.abs))
      private$leg.times.rel <- matrix(nrow = 1, ncol = ncol(private$leg.times.abs))
      for (i in 1:ncol(private$leg.times.abs)) {
        private$leg.times.rel[1, i] <- private$get.time.diff(time.vct[i], time.vct[i+1])
      }
      colnames(private$leg.times.rel) <- colnames(private$leg.times.abs)
    }
  )
)



#### Implementations of test of Bei (2024) ####

#' @title Perform the test of Bei (2024) for a given point
#'
#' @description This function performs the unconditional moment restriction test
#' as described in Bei (2024).
#'
#' @param r Result of the projection for which the test should be carried out.
#' @param c The projection matrix. For now, c is restricted to being an
#' elementary vector, i.e. c = (0, ...,0, 1, 0, ..., 0).
#' @param t The time point at which to evaluate theta.
#' @param par.space Matrix containing 2 columns and \eqn{d_\theta} rows, where
#' \eqn{d_\theta} is the dimension of the parameter space. The first column
#' represents the lower left corner of the parameter space, the second column
#' represents the upper right corner. At least for the time being, only
#' rectangular parameter spaces are allowed.
#' @param data Data frame on which to base the test.
#' @param hp List of hyperparameters needed.
#' @param verbose Boolean variable indicating whether to print updates of the
#' estimation process to the console.
#' @param inst.func.evals Matrix of precomputed instrumental function
#' evaluations for each observation in the data set. Used to speed up the
#' simulations. If \code{NULL}, the evaluations will be computed during
#' execution of this function. Default is \code{inst.func.evals = NULL}.
#' @param alpha The significance level at which to perform the test. Default is
#' \code{alpha = 0.95}.
#' @param parallel Flag for whether or not parallel computing should be used.
#' Default is \code{parallel = FALSE}.
#'
#' @importFrom foreach foreach
#'
#' @references Bei, X. (2024). Local linearization based subvector inference in
#' moment inequality models. Journal of Econometrics, 238(1),
#' 105549-. https://doi.org/10.1016/j.jeconom.2023.10554
#'
test.point_Bei <- function(r, c, t, par.space, data, hp, verbose = FALSE,
                           inst.func.evals = NULL, alpha = 0.95,
                           parallel = FALSE) {

  #### 0. Set some useful variables ####

  K.bar <- hp[["K.bar"]]
  n.cov <- sum(grepl("X[[:digit:]]+", colnames(data))) - 1
  n.param <- n.cov + 1
  delta.n <- hp[["delta.n"]]

  #### 1. Estimate the test statistic ####

  if (verbose) {
    message("\t Estimating test statistic...")
  }

  # Define initial value
  beta.init <- rowMeans(matrix(par.space[which(c == 0),], nrow = sum(c == 0)))

  # Get the test statistic
  out <- get.test.statistic(beta.init, data, par.space, t, hp, c, r,
                            inst.func.evals)
  Jnrh <- out[[1]]
  beta.hat <- out[[2]]

  # Define BetaI.r
  BetaI.r <- as.matrix(beta.hat)

  #### 2. Calculate the critical value for rn ####

  if (verbose) {
    message("\t Computing the initial critical value...")
  }

  cvLLn <- get.cvLLn(BetaI.r, data, t, hp, c, r, par.space, inst.func.evals,
                     alpha)

  #### 3. Obtain K.bar test statistics, starting from different values ####

  if (verbose & (Jnrh > cvLLn)) {
    message("\t Refinement step for critical value is skipped.")
  }

  if (Jnrh <= cvLLn) {
    if (verbose) {
      message(sprintf("\t Computing the K.bar (= %s) test statistics...", K.bar))
    }

    # Determine a set of K.bar initial values
    initial.values <- matrix(r, nrow = n.param, ncol = K.bar)
    for (idx in which(c == 0)) {
      initial.values[idx,] <- runif(K.bar, min = par.space[idx, 1], max = par.space[idx, 2])
    }

    # Initialize an object that will store all test statistic evaluations
    t.stat.evals <- matrix(nrow = K.bar, ncol = n.param + 1)
    colnames(t.stat.evals) <- c(paste0("X", 0:n.cov), "val")

    if (parallel) { # Using parallel computing

      t.stat.evals[] <-
        foreach(col.idx = 1:ncol(initial.values), .combine = 'rbind') %dopar% {

          source("lowLevelFunctions.R")

          # Select the parameter corresponding to this iteration
          beta.start <- initial.values[, col.idx]

          # Compute the test statistic
          out <- get.test.statistic(beta.start, data, par.space, t, hp, c, r,
                                    inst.func.evals)

          # Store the results
          c(out[[2]], out[[1]])
        }

    } else { # Using sequential computing

      # Compute \tilde{Beta}_I(r)
      for (col.idx in 1:ncol(initial.values)) {

        # Select the parameter corresponding to this iteration
        beta.start <- initial.values[, col.idx]

        # Compute the test statistic
        out <- get.test.statistic(beta.start, data, par.space, t, hp, c, r,
                                  inst.func.evals)

        # Store the results
        t.stat.evals[col.idx, ] <- c(out[[2]], out[[1]])
      }

    }

    # Reduce it to only the unique parameter vectors
    t.stat.evals <- t.stat.evals[which.unique(t.stat.evals[, 1:n.param]), , drop = FALSE]

    # Compute \hat{\beta}_I(r)
    t.stat.evals <- t.stat.evals[which(t.stat.evals[, n.param + 1] <= Jnrh + delta.n), , drop = FALSE]
  }

  #### 4. For each, calculate the critical value ####

  if (Jnrh <= cvLLn) {

    # If in the previous step no eligible values were found, this step is
    # skipped.
    if (nrow(t.stat.evals) == 0) {
      if (verbose) {
        message("\t No eligible minimizers found in refinement step")
      }

      # If eligible values were found, their critical value is computed.
    } else {
      if (verbose) {
        message("\t Recomputing the critical value...")
      }

      # Obtain the matrix of covariates in the correct format
      BetaI.r <- t(matrix(t.stat.evals[, 1:n.param], ncol = n.param))

      # Compute the critical value
      cvLLn <- get.cvLLn(BetaI.r, data, t, hp, c, r, par.space, inst.func.evals,
                         alpha)
    }
  }

  #### 5. Return rh, Jnrh, and cvLLn ####

  if (verbose) {
    message("\t Returning the results...")
  }

  # Note: 'as.numeric(cvLLn)' gets rid of the name of cvLLn, as it of type
  #       'named num'.
  c("theta" = r, "t.stat" = Jnrh, "crit.val" = as.numeric(cvLLn))
}

#' @title Perform the test of Bei (2024) simultaneously for multiple time
#' points.
#'
#' @description This function performs the unconditional moment restriction test
#' as described in Bei (2024). This function directly extends
#' \code{test.point_Bei} by allowing for pairs of moment restrictions over a
#' grid of time points.
#'
#' @param r Result of the projection for which the test should be carried out.
#' @param c The projection matrix. For now, c is restricted to being an
#' elementary vector, i.e. c = (0, ...,0, 1, 0, ..., 0).
#' @param t The time point at which to evaluate theta. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param par.space Matrix containing 2 columns and \eqn{d_\theta} rows, where
#' \eqn{d_\theta} is the dimension of the parameter space. The first column
#' represents the lower left corner of the parameter space, the second column
#' represents the upper right corner. At least for the time being, only
#' rectangular parameter spaces are allowed.
#' @param data Data frame on which to base the test.
#' @param hp List of hyperparameters needed.
#' @param verbose Boolean variable indicating whether to print updates of the
#' estimation process to the console.
#' @param inst.func.evals Matrix of precomputed instrumental function
#' evaluations for each observation in the data set. Used to speed up the
#' simulations. If \code{NULL}, the evaluations will be computed during
#' execution of this function. Default is \code{inst.func.evals = NULL}.
#' @param alpha The significance level at which to perform the test. Default is
#' \code{alpha = 0.95}.
#' @param parallel Flag for whether or not parallel computing should be used.
#' Default is \code{parallel = FALSE}.
#'
#' @importFrom foreach foreach
#'
#' @references Bei, X. (2024). Local linearization based subvector inference in
#' moment inequality models. Journal of Econometrics, 238(1),
#' 105549-. https://doi.org/10.1016/j.jeconom.2023.10554
#'
test.point_Bei_MT <- function(r, c, t, par.space, data, hp, verbose = FALSE,
                              inst.func.evals = NULL, alpha = 0.95,
                              parallel = FALSE) {

  #### 0. Set some useful variables ####

  # Hyperparameters
  K.bar <- hp[["K.bar"]]
  n.cov <- sum(grepl("X[[:digit:]]+", colnames(data))) - 1
  n.param <- n.cov + length(t)
  delta.n <- hp[["delta.n"]]

  #### 1. Estimate the test statistic ####

  if (verbose) {
    message("\t Estimating test statistic...")
  }

  # Define initial value
  beta.init <- rowMeans(matrix(par.space[which(c == 0),], nrow = sum(c == 0)))
  if (length(t) > 1) {
    beta.init[2:length(t)] <- par.space[2:length(t), 1] + 1
  }

  # Get the test statistic
  out <- get.test.statistic(beta.init, data, par.space, t, hp, c, r,
                            inst.func.evals)
  Jnrh <- out[[1]]
  beta.hat <- out[[2]]

  # Define BetaI.r
  BetaI.r <- as.matrix(beta.hat)

  #### 2. Calculate the critical value for rn ####

  if (verbose) {
    message("\t Computing the initial critical value...")
  }

  cvLLn <- get.cvLLn(BetaI.r, data, t, hp, c, r, par.space, inst.func.evals,
                     alpha)

  #### 3. Obtain K.bar test statistics, starting from different values ####

  if (verbose & (Jnrh > cvLLn)) {
    message("\t Refinement step for critical value is skipped.")
  }

  if (Jnrh <= cvLLn) {
    if (verbose) {
      message(sprintf("\t Computing the K.bar (= %s) test statistics...", K.bar))
    }

    # Determine a set of K.bar initial values
    initial.values <- matrix(r, nrow = n.param, ncol = K.bar)
    for (idx in which(c == 0)) {
      initial.values[idx,] <- runif(K.bar, min = par.space[idx, 1], max = par.space[idx, 2])
    }

    # Initialize an object that will store all test statistic evaluations
    t.stat.evals <- matrix(nrow = K.bar, ncol = n.param + 1)
    if (length(t) == 1) {
      colnames(t.stat.evals) <- c(paste0("X", 0:n.cov), "val")
    } else {
      colnames(t.stat.evals) <- c(paste0("X0.", t),
                                  paste0("X", 1:n.cov),
                                  "val")
    }

    if (parallel) { # Using parallel computing

      t.stat.evals[] <-
        foreach(col.idx = 1:ncol(initial.values), .combine = 'rbind') %dopar% {

          source("lowLevelFunctions.R")

          # Select the parameter corresponding to this iteration
          beta.start <- initial.values[, col.idx]

          # Compute the test statistic
          out <- get.test.statistic(beta.start, data, par.space, t, hp, c, r,
                                    inst.func.evals)

          # Store the results
          c(out[[2]], out[[1]])
        }

    } else { # Using sequential computing

      # Compute \tilde{Beta}_I(r)
      for (col.idx in 1:ncol(initial.values)) {

        # Select the parameter corresponding to this iteration
        beta.start <- initial.values[, col.idx]

        # Compute the test statistic
        out <- get.test.statistic(beta.start, data, par.space, t, hp, c, r,
                                  inst.func.evals)

        # Store the results
        t.stat.evals[col.idx, ] <- c(out[[2]], out[[1]])
      }

    }

    # Reduce it to only the unique parameter vectors
    t.stat.evals <- t.stat.evals[which.unique(t.stat.evals[, 1:n.param]), , drop = FALSE]

    # Compute \hat{\beta}_I(r)
    t.stat.evals <- t.stat.evals[which(t.stat.evals[, n.param + 1] <= Jnrh + delta.n), , drop = FALSE]
  }

  #### 4. For each, calculate the critical value ####

  if (Jnrh <= cvLLn) {

    # If in the previous step no eligible values were found, this step is
    # skipped.
    if (nrow(t.stat.evals) == 0) {
      if (verbose) {
        message("\t No eligible minimizers found in refinement step")
      }

      # If eligible values were found, their critical value is computed.
    } else {
      if (verbose) {
        message("\t Recomputing the critical value...")
      }

      # Obtain the matrix of covariates in the correct format
      BetaI.r <- t(matrix(t.stat.evals[, 1:n.param], ncol = n.param))

      # Compute the critical value
      cvLLn <- get.cvLLn(BetaI.r, data, t, hp, c, r, par.space, inst.func.evals,
                         alpha)
    }
  }

  #### 5. Return rh, Jnrh, and cvLLn ####

  if (verbose) {
    message("\t Returning the results...")
  }

  # Note: 'as.numeric(cvLLn)' gets rid of the name of cvLLn, as it of type
  #       'named num'.
  c("theta" = r, "t.stat" = Jnrh, "crit.val" = as.numeric(cvLLn))
}


#### Low level functions: general utility functions ####

#' @title Tests whether two parameter vectors are approximately equal
#'
#' @description
#' Function to check whether the given parameters are approximately equal, where
#' approximate equality is defined as the norm of the differences being smaller
#' than the given tolerance.
#'
#' @param par1 First parameter
#' @param par2 Second parameter
#' @param tol Tolerance. Differences between \code{par1} and \code{par2} smaller
#' than \code{tol} will be neglected. Default value is \code{tol = 10^(-4)}.
#'
#' @noRd
#'
par.equal <- function(par1, par2, tol = 10^(-4)) {
  distance <- sqrt(sum((par1 - par2)^2))
  distance < tol
}

#' @title Returns the indices of the (approximately) unique parameter vectors
#'
#' @description
#' Returns the indices of the approximately unique parameter vectors in a
#' matrix, where unicity is defined as there not existing another parameter
#' vector that is equal up until a difference in norm less than the given
#' tolerance.
#'
#' @param par.mat Matrix of parameter vectors.
#' @param tol Tolerance. Default value is \code{tol = 10^(-4)}.
#'
#' @noRd
#'
which.unique <- function(par.mat, tol = 10^(-4)) {

  # If par.mat is a vector, return 1
  if (!is.matrix(par.mat)) {
    return(c(1))
  }

  # Initialize a vector that will store all unique column indices
  idx.unique <- c()

  # For all but the last row, check whether there are duplicate rows that
  # follow.
  if (nrow(par.mat) > 1) {
    for (i in 1:(nrow(par.mat) - 1)) {
      dupe <- FALSE
      for (j in (i+1):nrow(par.mat)) {
        if (par.equal(par.mat[i,], par.mat[j,], tol = tol)) {
          dupe <- TRUE
        }
      }

      if (!dupe) {
        idx.unique <- c(idx.unique, i)
      }
    }
  }

  # The last element is always unique
  idx.unique <- c(idx.unique, nrow(par.mat))

  # Return the result
  idx.unique
}

#' @title Checks if a directory exists and, if necessary, creates it.
#'
#' @description
#' Function that checks whether a given directory exists and, if necessary,
#' creates it.
#'
#' @param dir.name Name of the directory whose existence is to be checked.
#' @param path Path to the parent directory in which the given directory should
#' exist. Default value is \code{path = NULL}, in which case the working
#' directory is used.
#'
#' @noRd
#'
check_create.dir <- function(dir.name, path = NULL) {

  # Path to directory
  if (!is.null(path)) {
    dir.path <- paste(c(path, dir.name), collapse = "/")
  } else {
    dir.path <- dir.name
  }

  # If the directory does not exist yet, create it
  if (!dir.exists(dir.path)) {
    dir.create(dir.path)
  }
}

#### Low level functions: data generation function ####

#' @title Test whether covariates lie inside {1}^d.d x [0.5, 1.5]^d.c
#'
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {1}^d.d x
#' [0.5, 1.5]^d.c. This function is used in generating data according to some
#' DGPs in the function 'generateData.R'.
#'
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#'
#' @noRd
#'
inRegion1 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  check.in.region <- function(row) {
    all(0.5 <= row[type.cov == "c"]) & all(row[type.cov == "c"] <= 1.5) &
      all(row[type.cov == "d"] == 1)
  }
  apply(X[, -1, drop = FALSE], 1, check.in.region)
}

#' @title Test whether covariates lie inside {1}^d.d x [-Inf, -1]^d
#'
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {1}^d.d x
#' [-Inf, 1]^d. This function is used in generating data according to some DGPs
#' in the function 'generateData.R'.
#'
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#'
#' @noRd
#'
inRegion2 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  check.in.region <- function(row) {
    all(row[type.cov == "c"] <= -1) & all(row[type.cov == "d"] == 1)
  }
  apply(X[, -1, drop = FALSE], 1, check.in.region)
}

#' @title Test whether covariates lie inside {0}^d.d x ([-1, 0] x [0, 1])^(d/2)
#'
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {0}^d.d x
#' ([-1, 0] x [0, 1])^(d/2). This function is used in generating data according
#' to some DGPs in the function 'generateData.R'. This function is also used in
#' the implementation of \code{'inRegion4.R'}.
#'
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#'
#' @noRd
#'
inRegion3 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  check.in.region <- function(row) {
    lbs <- rep(c(-1, 0), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    ubs <- rep(c(0, 1), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    equal.d <- all(row[type.cov != "c"] == 1)
    equal.c <- all(lbs <= row[type.cov == "c"]) & all(row[type.cov == "c"] <= ubs)
    equal.d & equal.c
  }
  apply(X[, -1, drop = FALSE], 1, check.in.region)
}

#' @title Test whether covariates lie inside {0}^d.d x ([0, 1] x [-1, 0])^(d/2)
#'
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {0}^d.d x
#' ([0, 1] x [-1, 0])^(d/2). This function is used in generating data according
#' to some DGPs in the function 'generateData.R'.
#'
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#'
#' @noRd
#'
inRegion4 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  inRegion3(-X, type.cov)
}

#' @title Test whether covariates lie inside {0}^d.d x ([-0.5, 0.5] x
#' [-2, 2])^(d/2)
#'
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {0}^d.d x
#' ([-0.5, 0.5] x [-2, 2])^(d/2). This function is used in generating data
#' according to some DGPs in the function 'generateData_add.R'.
#'
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#'
#' @noRd
#'
inRegion5 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  check.in.region <- function(row) {
    lbs <- rep(c(-0.5, -2), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    ubs <- rep(c(0.5, 2), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    equal.d <- all(row[type.cov != "c"] == 1)
    equal.c <- all(lbs <= row[type.cov == "c"]) & all(row[type.cov == "c"] <= ubs)
    equal.d & equal.c
  }
  apply(X[, -1, drop = FALSE], 1, check.in.region)
}

#' @title Test whether covariates lie inside {0}^d.d x ([-2, 2] x
#' [-0.5, 0.5])^(d/2)
#'
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {0}^d.d x
#' ([-2, 2] x [-0.5, 0.5])^(d/2). This function is used in generating data
#' according to some DGPs in the function 'generateData_add.R'.
#'
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#'
#' @noRd
#'
inRegion6 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  check.in.region <- function(row) {
    lbs <- rep(c(-2, -0.5), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    ubs <- rep(c(2, 0.5), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    equal.d <- all(row[type.cov != "c"] == 1)
    equal.c <- all(lbs <= row[type.cov == "c"]) & all(row[type.cov == "c"] <= ubs)
    equal.d & equal.c
  }
  apply(X[, -1, drop = FALSE], 1, check.in.region)
}

#' @title Generates a data set according to the specified arguments.
#'
#' @description This function generates a data set according to the specified
#' arguments.
#'
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options List of additional arguments.
#' @param plot_data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot_data = FALSE}.
#'
#' @import stats
#' @importFrom graphics hist
#' @importFrom copula frankCopula rCopula
#'
#' @noRd
#'
generateData <- function(beta.true, n, n.cov, options, plot_data = FALSE) {

  # Extract the necessary hyperparameters
  if (options[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options[["DGP"]]

  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R',
  #               'simFuncWrapper.R' and 'simulate1D.CCDF.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }

  # For the data generation, we can make the following shortcut.
  beta.true <- beta.true(0)

  # Subset beta.true to the correct number of parameters
  beta.true <- beta.true[1:(n.cov + 1)]

  # Generate the intercept and the covariates
  X <- rep(1, n)
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, rnorm(n))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, rbinom(n, 1, 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)

  # Determine which covariate vectors lie in the predefined regions.
  idxs.region1 <- which(inRegion1(X, type.cov))
  idxs.region2 <- which(inRegion2(X, type.cov))
  idxs.region3 <- which(inRegion3(X, type.cov))
  idxs.region4 <- which(inRegion4(X, type.cov))

  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # Independence, C ~ Unif

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    C <- pmax(runif(n, min(T), max(T)), runif(n, min(T), max(T)))

  } else if (DGP %% 20 == 2) { # Positive dependence, C ~ Unif

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Define quantile function of C
    FC_given_X_inv <- function(u2) {u2 * (max(T) - min(T)) + min(T) + 2}
    C <- FC_given_X_inv(u2)

  } else if (DGP %% 20 == 3) { # Negative dependence, C ~ Unif

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Define quantile function of C
    FC_given_X_inv <- function(u2) {u2 * (max(T) - min(T)) + min(T) + 4}
    C <- FC_given_X_inv(u2)

  } else if (DGP %% 20 == 4) { # Independence, C ~ AFT_ll

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- gamma[1] - 1.3
    C <- Lambda_inverse_AFT_ll(runif(n)) - X %*% gamma

  } else if (DGP %% 20 == 5) { # Positive dependence, C ~ AFT_ll

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- gamma[1] - 1.3
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma

  } else if (DGP %% 20 == 6) { # Negative dependence, C ~ AFT_ll

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- gamma[1] - 1.3
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma

  } else if (DGP %% 20 == 7) { # Independence, C ~ AFT_ll, high cens

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    C <- Lambda_inverse_AFT_ll(runif(n)) - X %*% gamma

  } else if (DGP %% 20 == 8) { # Positive dependence, C ~ AFT_ll, high cens

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma

  } else if (DGP %% 20 == 9) { # Negative dependence, C ~ AFT_ll, high cens

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma

  } else if (DGP %% 20 == 10) { # Independence, C ~ AFT_ll, high cens, regions with less cens

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region2 <- beta.true
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region2[1] <- beta.true[1] - 2

    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion2(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(runif(length(idxs.region1))) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region2] <- Lambda_inverse_AFT_ll(runif(length(idxs.region2))) - X[idxs.region2, ] %*% gamma.region2

  } else if (DGP %% 20 == 11) { # Pos. dep., C ~ AFT_ll, high cens, regions with less cens

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region2 <- beta.true
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region2[1] <- beta.true[1] - 2

    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion2(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region2] <- Lambda_inverse_AFT_ll(u2[idxs.region2]) - X[idxs.region2, ] %*% gamma.region2

  } else if (DGP %% 20 == 12) { # Neg. dep., C ~ AFT_ll, high cens, regions with less cens

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region2 <- beta.true
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region2[1] <- beta.true[1] - 2

    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion2(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region2] <- Lambda_inverse_AFT_ll(u2[idxs.region2]) - X[idxs.region2, ] %*% gamma.region2

  } else if (DGP %% 20 == 13) { # Positive dependence, cens comparable to DGP = 11

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- beta.true[1] - 0.115
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma

  } else if (DGP %% 20 == 14) { # Negative dependence, cens comparable to DGP = 12

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- beta.true[1] - 0.2
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma

  } else if (DGP %% 20 == 15) { # Independence, C ~ AFT_ll, high cens, regions with less cens 2

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.43
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region4[1] <- beta.true[1] - 2

    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(runif(length(idxs.region1))) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(runif(length(idxs.region4))) - X[idxs.region4, ] %*% gamma.region4

  } else if (DGP %% 20 == 16) { # Pos. dep., C ~ AFT_ll, high cens, regions with less cens 2

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.15
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region4[1] <- beta.true[1] - 2

    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4

  } else if (DGP %% 20 == 17) { # Neg. dep., C ~ AFT_ll, high cens, regions with less cens 2

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region4[1] <- beta.true[1] - 2

    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4

  } else if (DGP %% 20 == 18) { # Independence, C ~ AFT_ll, high cens, regions with less cens 3

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region3 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.58
    gamma.region3[1] <- beta.true[1] - 3
    gamma.region4[1] <- beta.true[1] - 2

    idxs.outside <- which(!inRegion3(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region3] <- Lambda_inverse_AFT_ll(runif(length(idxs.region3))) - X[idxs.region3, ] %*% gamma.region3
    C[idxs.region4] <- Lambda_inverse_AFT_ll(runif(length(idxs.region4))) - X[idxs.region4, ] %*% gamma.region4

  } else if (DGP %% 20 == 19) { # Pos. dep., C ~ AFT_ll, high cens, regions with less cens 3

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region3 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.25
    gamma.region3[1] <- beta.true[1] - 3
    gamma.region4[1] <- beta.true[1] - 2

    idxs.outside <- which(!inRegion3(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region3] <- Lambda_inverse_AFT_ll(u2[idxs.region3]) - X[idxs.region3, ] %*% gamma.region3
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4

  }

  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)

  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))

  # Histogram of the observed times
  if (plot_data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }

  # Return the results
  data
}

#' @title Additional data generating function.
#'
#' @description This function generates a data set according to the specified
#' arguments, like 'generateData.R' above. It differs from the aforementioned
#' function in that some DGP's are slightly different
#'
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options List of additional arguments.
#' @param plot_data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot_data = FALSE}.
#'
#' @import stats
#' @importFrom graphics hist
#' @importFrom copula frankCopula rCopula
#'
#' @noRd
#'
generateData_add <- function(beta.true, n, n.cov, options, plot_data = FALSE) {

  # Extract the necessary hyperparameters
  if (options[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options[["DGP"]]

  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R' and
  #               'simFuncWrapper.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }

  # Because, at least for now, beta.true only depends on t in its first element
  # (corresponding to the intercept) and moreover, this dependence is of the
  # form beta_0(t) = t + a, we can make the following short-cut.
  beta.true <- beta.true(0)

  # Subset beta.true to the correct number of parameters
  beta.true <- beta.true[1:(n.cov + 1)]

  # Generate the intercept and the covariates
  X <- rep(1, n)
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, runif(n, -3, 3))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, rbinom(n, 1, 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)

  # Determine which covariate vectors lie in the predefined regions.
  idxs.region1 <- which(inRegion1(X, type.cov))
  idxs.region2 <- which(inRegion2(X, type.cov))
  idxs.region3 <- which(inRegion3(X, type.cov))
  idxs.region4 <- which(inRegion4(X, type.cov))
  idxs.region5 <- which(inRegion5(X, type.cov))
  idxs.region6 <- which(inRegion6(X, type.cov))

  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # Independence, region 5

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- beta.true
    gamma[1] <- beta.true[1] + 0.4
    gamma.region5[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion5(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(runif(length(idxs.region5))) - X[idxs.region5, ] %*% gamma.region5

  } else if (DGP %% 20 == 2) { # Positive dependence, region 5

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- beta.true
    gamma[1] <- beta.true[1] + 0
    gamma.region5[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion5(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(u2[idxs.region5]) - X[idxs.region5, ] %*% gamma.region5

  } else if (DGP %% 20 == 3) { # Negative dependence, region 5

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region5[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion5(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(u2[idxs.region5]) - X[idxs.region5, ] %*% gamma.region5

  } else if (DGP %% 20 == 4) { # Independence, region 6

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.4
    gamma.region6[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region6] <- Lambda_inverse_AFT_ll(runif(length(idxs.region6))) - X[idxs.region6, ] %*% gamma.region6

  } else if (DGP %% 20 == 5) { # Positive dependence, region 6

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.05
    gamma.region6[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region6] <- Lambda_inverse_AFT_ll(u2[idxs.region6]) - X[idxs.region6, ] %*% gamma.region6

  } else if (DGP %% 20 == 6) { # Negative dependence, region 6

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region6[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region6] <- Lambda_inverse_AFT_ll(u2[idxs.region6]) - X[idxs.region6, ] %*% gamma.region6

  } else if (DGP %% 20 == 7) { # Independence, regions 5 and 6

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.8
    gamma.region5[1] <- gamma.region6[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion5(X, type.cov) & !inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(runif(length(idxs.region5))) - X[idxs.region5, ] %*% gamma.region5
    C[idxs.region6] <- Lambda_inverse_AFT_ll(runif(length(idxs.region6))) - X[idxs.region6, ] %*% gamma.region6

  } else if (DGP %% 20 == 8) { # Positive dependence, regions 5 and 6

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region5[1] <- gamma.region6[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion5(X, type.cov) & !inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(u2[idxs.region5]) - X[idxs.region5, ] %*% gamma.region5
    C[idxs.region6] <- Lambda_inverse_AFT_ll(u2[idxs.region6]) - X[idxs.region6, ] %*% gamma.region6

  } else if (DGP %% 20 == 9) { # Negative dependence, regions 5 and 6

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.7
    gamma.region5[1] <- gamma.region6[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion5(X, type.cov) & !inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(u2[idxs.region5]) - X[idxs.region5, ] %*% gamma.region5
    C[idxs.region6] <- Lambda_inverse_AFT_ll(u2[idxs.region6]) - X[idxs.region6, ] %*% gamma.region6

  } else if (DGP %% 20 == 10) { # Independence, regions 1 + 4

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region1[1] <- beta.true[1] - 5
    gamma.region4[1] <- beta.true[1] - 6

    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(runif(length(idxs.region1))) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(runif(length(idxs.region4))) - X[idxs.region4, ] %*% gamma.region4

  } else if (DGP %% 20 == 11) { # Positive dependence, regions 1 + 4

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1]
    gamma.region1[1] <- beta.true[1] - 4
    gamma.region4[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4

  } else if (DGP %% 20 == 12) { # Negative dependence, regions 1 + 4

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.1
    gamma.region1[1] <- beta.true[1] - 6
    gamma.region4[1] <- beta.true[1] - 7

    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4

  } else if (DGP %% 20 == 13) { # Independence, regions 3 + 4

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region3 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region3[1] <- beta.true[1] - 5
    gamma.region4[1] <- beta.true[1] - 6

    idxs.outside <- which(!inRegion3(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region3] <- Lambda_inverse_AFT_ll(runif(length(idxs.region3))) - X[idxs.region3, ] %*% gamma.region3
    C[idxs.region4] <- Lambda_inverse_AFT_ll(runif(length(idxs.region4))) - X[idxs.region4, ] %*% gamma.region4

  } else if (DGP %% 20 == 14) { # Positive dependence, regions 3 + 4

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region3 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1]
    gamma.region3[1] <- beta.true[1] - 5
    gamma.region4[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion3(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region3] <- Lambda_inverse_AFT_ll(u2[idxs.region3]) - X[idxs.region3, ] %*% gamma.region3
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4

  } else if (DGP %% 20 == 15) { # Negative dependence, regions 3 + 4

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma.region3 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1]
    gamma.region3[1] <- beta.true[1] - 5
    gamma.region4[1] <- beta.true[1] - 5

    idxs.outside <- which(!inRegion3(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region3] <- Lambda_inverse_AFT_ll(u2[idxs.region3]) - X[idxs.region3, ] %*% gamma.region3
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4

  } else if (DGP %% 20 == 16) { # Independence, C ~ AFT_ll, high cens. ~ DGP 7 in generateData

    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    C <- Lambda_inverse_AFT_ll(runif(n)) - X %*% gamma

  } else if (DGP %% 20 == 17) { # Pos. dep., C ~ AFT_ll, high cens. ~ DGP 13 in generateData

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- beta.true[1] - 0.115
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma

  } else if (DGP %% 20 == 18) { # Neg. dep., C ~ AFT_ll, high cens. ~ DGP 14 in generateData

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true

    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- beta.true[1] - 0.2
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma

  }

  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)

  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))

  # Histogram of the observed times
  if (plot_data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }

  # Return the results
  data
}

#' @title Data generation function for the main simulation.
#'
#' @description This function generates a data set according to the specified
#' arguments.
#'
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options List of additional arguments.
#' @param H0.inv Inverse of the intercept (function of time).
#' @param plot_data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot_data = FALSE}.
#'
#' @import stats
#' @importFrom graphics hist
#' @importFrom copula frankCopula rCopula
#'
#' @noRd
#'
generateData_simMain <- function(beta.true, n, n.cov, options, H0.inv,
                                 plot_data = FALSE) {

  # Extract the necessary hyperparameters
  if (options[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options[["DGP"]]

  # Subset the parameter vector to the first n.cov covariate effects.
  if (is(beta.true, "function")) {
    beta.true <- beta.true(0)[2:(n.cov + 1)]
  }

  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R',
  #               'simFuncWrapper.R' and 'simulate1D.CCDF.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }

  # Generate the intercept and the covariates
  X <- rep(1, n)
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, rnorm(n))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, rbinom(n, 1, 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)

  # Get matrix of just the covariates (no intercept)
  X.noint <- X[, -1]

  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # For AFT_ll: Independence, C ~ Exp, ~25% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.1
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 2) { # For AFT_ll: Pos. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.15
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 3) { # For AFT_ll: Neg. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.07
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 4) { # For AFT_ll: Independence, C ~ Exp, ~65% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.8
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 5) { # For AFT_ll: Pos. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.65
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 6) { # For AFT_ll: Neg. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 1.2
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 7) { # For Cox_wb: Independence, C ~ Exp, ~25% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.22
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 8) { # For Cox_wb: Pos. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.3
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 9) { # For Cox_wb: Neg. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.17
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 10) { # For Cox_wb: Independence, C ~ Exp, ~65% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 1.3
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 11) { # For Cox_wb: Pos. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 1
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 12) { # For Cox_wb: Neg. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 1.6
    C <- (-1/lambda.c) * log(1 - u2)

  }  else if (DGP %% 20 == 13) { # For AFT_ll: Independence, C ~ Exp, ~2% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.002
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 14) { # For AFT_ll: Pos. dep., C ~ Exp, ~2% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.006
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 15) { # For AFT_ll: Neg. dep., C ~ Exp, ~2% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.0005
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 16) { # For Cox_wb: Independence, C ~ Exp, ~2% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.005
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 17) { # For Cox_wb: Pos. dep., C ~ Exp, ~2% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.03
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 18) { # For Cox_wb: Neg. dep., C ~ Exp, ~2% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.003
    C <- (-1/lambda.c) * log(1 - u2)

  }

  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)

  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))

  # Histogram of the observed times
  if (plot_data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }

  # Return the results
  data
}

#' @title Data generation function for the additional simulation.
#'
#' @description This function generates a data set according to the specified
#' arguments. As opposed to \code{generateData_simMain.R}, it generates the
#' covariates in a way that they are depedendent. (Achieved through the use of
#' a copula).
#'
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options List of additional arguments.
#' @param H0.inv Inverse of the intercept (function of time).
#' @param plot_data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot_data = FALSE}.
#'
#' @import stats
#' @importFrom graphics hist
#' @importFrom copula frankCopula normalCopula rCopula
#'
#' @noRd
#'
generateData_simAdd <- function(beta.true, n, n.cov, options, H0.inv,
                                plot_data = FALSE) {

  # Extract the necessary hyperparameters
  if (options[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options[["DGP"]]

  # Subset the parameter vector to the first n.cov covariate effects.
  if (is(beta.true, "function")) {
    beta.true <- beta.true(0)[2:(n.cov + 1)]
  }

  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R',
  #               'simFuncWrapper.R' and 'simulate1D.CCDF.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }

  # Throw an error if specified number of covariates is not equal to two. These
  # cases are not implemented in this function.
  if (n.cov != 2) {
    stop("The specified number of covariates should equal 2.")
  }

  # Generate the intercept and the covariates
  X <- rep(1, n)
  X.U <- rCopula(n, normalCopula(0.8))
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, qnorm(X.U[,i]))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, as.numeric(X.U[,i] >= 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)

  # Get matrix of just the covariates (no intercept)
  X.noint <- X[, -1]

  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # For AFT_ll: Independence, C ~ Exp, ~30% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.13
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 2) { # For AFT_ll: Pos. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.17
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 3) { # For AFT_ll: Neg. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.10
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 4) { # For AFT_ll: Independence, C ~ Exp, ~65% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 1
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 5) { # For AFT_ll: Pos. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.6
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 6) { # For AFT_ll: Neg. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 1.3
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 7) { # For Cox_wb: Independence, C ~ Exp, ~25% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.25
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 8) { # For Cox_wb: Pos. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.31
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 9) { # For Cox_wb: Neg. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.19
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 10) { # For Cox_wb: Independence, C ~ Exp, ~65% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 1.3
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 11) { # For Cox_wb: Pos. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.9
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 12) { # For Cox_wb: Neg. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 1.5
    C <- (-1/lambda.c) * log(1 - u2)

  }  else if (DGP %% 20 == 13) { # For AFT_ll: Independence, C ~ Exp, ~2% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.002
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 14) { # For AFT_ll: Pos. dep., C ~ Exp, ~2% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.01
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 15) { # For AFT_ll: Neg. dep., C ~ Exp, ~2% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.04
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 16) { # For Cox_wb: Independence, C ~ Exp, ~2% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.01
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 17) { # For Cox_wb: Pos. dep., C ~ Exp, ~2% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.03
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 18) { # For Cox_wb: Neg. dep., C ~ Exp, ~2% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.005
    C <- (-1/lambda.c) * log(1 - u2)

  }

  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)

  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))

  # Histogram of the observed times
  if (plot_data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }

  # Return the results
  data
}

#' @title Data generation function for the simulation regarding misspecification.
#'
#' @description This function generates a data set according to the specified
#' arguments. It is mostly a copy-paste from the function generateData_simMain.R
#' above, but the censoring distribution slightly addapted in order to control
#' the percentage of censored observations.
#'
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options.data.gen List of additional arguments.
#' @param H0.inv Inverse of the intercept (function of time).
#' @param plot_data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot_data = FALSE}.
#'
#' @import stats
#' @importFrom graphics hist
#' @importFrom copula frankCopula rCopula
#'
#' @noRd
#'
generateData_simMiss <- function(beta.true, n, n.cov, options.data.gen, H0.inv,
                                 plot_data = FALSE) {

  # Parameter used for development of this function. Set to TRUE MANUALLY if
  # desired.
  # NOTE: When the function is run as is with 'test.mode = TRUE', it will fail!
  test.mode <- FALSE

  # Extract the necessary hyperparameters
  if (options.data.gen[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options.data.gen[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options.data.gen[["DGP"]]

  # Subset the parameter vector to the first n.cov covariate effects.
  if (is(beta.true, "function")) {
    beta.true <- beta.true(0)[2:(n.cov + 1)]
  }

  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R',
  #               'simFuncWrapper.R' and 'simulate1D.CCDF.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }

  # Generate the intercept and the covariates
  X <- rep(1, n)
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, rnorm(n))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, rbinom(n, 1, 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)

  # Get matrix of just the covariates (no intercept)
  X.noint <- X[, -1]

  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # For AFT_ll: Independence, C ~ Exp, ~25% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.12
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 2) { # For AFT_ll: Pos. dep., C ~ Unif, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.18
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 4) { # For AFT_ll: Independence, C ~ Exp, ~65% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 1
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 5) { # For AFT_ll: Pos. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.75
    C <- (-1/lambda.c) * log(1 - u2)

  }

  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)

  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))

  # Histogram of the observed times
  if (plot_data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }

  # Return the results
  data
}

#' @title Data generation function for the simulations under many covariates.
#'
#' @description This function generates a data set according to the specified
#' arguments.
#'
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options List of additional arguments.
#' @param H0.inv Inverse of the intercept (function of time).
#' @param plot_data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot_data = FALSE}.
#'
#' @import stats
#' @importFrom graphics hist
#' @importFrom copula frankCopula rCopula
#'
#' @noRd
#'
generateData_simManyCov <- function(beta.true, n, n.cov, options, H0.inv,
                                    plot_data = FALSE) {

  # Extract the necessary hyperparameters
  if (options[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options[["DGP"]]

  # Subset the parameter vector to the first n.cov covariate effects.
  if (is(beta.true, "function")) {
    beta.true <- beta.true(0)[2:(n.cov + 1)]
  }

  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R',
  #               'simFuncWrapper.R' and 'simulate1D.CCDF.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }

  # Generate the intercept and the covariates
  X <- rep(1, n)
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, rnorm(n))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, rbinom(n, 1, 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)

  # Get matrix of just the covariates (no intercept)
  X.noint <- X[, -1]

  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # For AFT_ll: Independence, C ~ Exp, ~25% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.35
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 2) { # For AFT_ll: Pos. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.45
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 3) { # For AFT_ll: Neg. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.25
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 4) { # For AFT_ll: Independence, C ~ Exp, ~65% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 5
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 5) { # For AFT_ll: Pos. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 3.5
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 6) { # For AFT_ll: Neg. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 6
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 7) { # For Cox_wb: Independence, C ~ Exp, ~25% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.7
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 8) { # For Cox_wb: Pos. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.9
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 9) { # For Cox_wb: Neg. dep., C ~ Exp, ~25% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 0.55
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 10) { # For Cox_wb: Independence, C ~ Exp, ~65% cens.

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 8
    C <- (-1/lambda.c) * log(1 - runif(n))

  } else if (DGP %% 20 == 11) { # For Cox_wb: Pos. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 6
    C <- (-1/lambda.c) * log(1 - u2)

  } else if (DGP %% 20 == 12) { # For Cox_wb: Neg. dep., C ~ Exp, ~65% cens.

    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)

    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]

    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))

    # Latent censoring time C
    lambda.c <- 9
    C <- (-1/lambda.c) * log(1 - u2)

  }

  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)

  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))

  # Histogram of the observed times
  if (plot_data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }

  # Return the results
  data
}

#### Low level functions: link functions ####

#' @title Link function (Cox model)
#'
#' @description
#' This function defines the Cox PH link function.
#'
#' @param t time parameter.
#'
Lambda_Cox_wb <- function(t) {
  1 - exp(-exp(t))
}

#' @title Derivative of link function (Cox model)
#'
#' @description
#' This function defines the derivative of the Cox PH link function.
#'
#' @param t time parameter.
#'
dLambda_Cox_wb <- function(t) {
  exp(t - exp(t))
}

#' @title Inverse of link function (Cox model)
#'
#' @description
#' This function defines the inverse of the Cox PH link function.
#'
#' @param p probability.
#'
Lambda_inverse_Cox_wb <- function(p) {
  log(-log(1-p))
}

#' @title Link function (AFT model)
#'
#' @description
#' This function defines the AFT link function.
#'
#' @param t time parameter.
#'
Lambda_AFT_ll <- function(t) {
  1 - 1/(1 + exp(t))
}

#' @title Derivative of link function (AFT model)
#'
#' @description
#' This function defines the derivative of the AFT link function.
#'
#' @param t time parameter.
#'
dLambda_AFT_ll <- function(t) {
  exp(t)/(1 + exp(t))^2
}

#' @title Inverse of link function (AFT model)
#'
#' @description
#' This function defines the inverse of the AFT link function.
#'
#' @param p probability.
#'
Lambda_inverse_AFT_ll <- function(p) {
  log(p / (1 - p))
}

#### Low level functions: family of instrumental functions ####

#' @title Normalize the covariates of a data set to lie in the unit interval by
#' scaling based on the ranges of the covariates.
#'
#' @description This function normalized the covariates in the data to lie in
#' the unit interval based on either the empirical or known ranges of the
#' covariates. It is useful to perform this step when defining the instrumental
#' functions later on. This function is used in \code{G.box.R}, \code{G.spline.R}
#' and by extension in \code{G.cd.R}.
#'
#' @param data (optional) Data set to be used to construct the normalizing
#' transformation. Default is \code{data = NULL}.
#' @param x (optional) Vector of covariates to be normalized alongside the data.
#' Default is \code{x = NULL}.
#' @param cov.ranges (optional) Matrix that specifies the range of each of the
#' covariates in the data set. Each column corresponds to a covariate. The first
#' row contains the lower bound, the second row contains the upper bound.
#' If not supplied, the data will be normalized based on the minimum and maximum
#' detected values. If supplied, the non data-dependent transformation function
#' listed in the appendix of Andrews, Shi 2013 will be used. Default is
#' \code{cov.ranges = NULL}.
#' @param norm.cov.out (optional) The output of a previous call to this function.
#' Can be used to speed up computation. If both \code{data} and
#' \code{norm.cov.out} are supplied to the function, this method will throw an
#' error. Default is \code{norm.cov.out = NULL}.
#' @param idxs.c (optional) Vector of indices of covariates that are continuous.
#' Note that that indices are relative to the covariate vector, not the full
#' data set. Default value is \code{idxs.c = "all"}, which indicates that all
#' elements should be regarded as continuous. If \code{idxs.c = NULL}, all
#' elements are regarded as discrete.
#' @param ... Allows easier interchangeability between covariate normalization
#' functions. All arguments specified under \code{...} will be ignored.
#'
#' @references Andrews, D.W.K. and Shi, X. (2013). Inference based on
#' confitional moment inequalities. Econometrica. 81(2):609-666.
#'
normalize.covariates <- function(data = NULL, x = NULL, cov.ranges = NULL,
                                 idxs.c = "all", norm.cov.out = NULL, ...) {

  # Precondition checks
  if (is.null(data) & is.null(norm.cov.out)) {
    stop("Either data or norm.cov.out should be supplied to this function.")
  }
  if (!is.null(data) & !is.null(norm.cov.out)) {
    stop("Ambiguous function arguments: both data and norm.cov.out are supplied.")
  }

  # Extract the covariates from the data set, if applicable. Else extract the
  # necessary parameters from the previous function call.
  if (!is.null(data)) {

    # Get all covariates names
    cov.idxs <- which(grepl("X[[:digit:]]+", colnames(data)))[-1]
    covariate.names <- colnames(data)[cov.idxs]

    # Only retain names of continuous covariates
    idxs.c <- if (all(idxs.c == "all")) {1:length(covariate.names)} else {idxs.c}
    cont.cov.names <- covariate.names[idxs.c]
  } else {
    covariate.names <- norm.cov.out$covariate.names
    cont.cov.names <- norm.cov.out$cont.cov.names
  }

  # If supplied, rename the entries of x.
  if (!is.null(x)) {
    names(x) <- covariate.names
  }

  # Initialize object that will store the data with normalized covariates
  normalized.data <- data
  x.c.norm <- NULL

  # Compute the minimum and maximum value of each continuous covariates based on
  # the data. If the output of a previous call to this function was provided,
  # this step can be skipped.
  X.ranges <- matrix(nrow = 2, ncol = length(cont.cov.names))
  colnames(X.ranges) <- cont.cov.names
  if (is.null(norm.cov.out)) {
    for (cov.name in cont.cov.names) {
      if (is.null(cov.ranges)) {
        X.ranges[1, cov.name] <- min(data[, cov.name])
        X.ranges[2, cov.name] <- max(data[, cov.name])
      } else {
        X.ranges[1, cov.name] <- cov.ranges[1, cov.name]
        X.ranges[2, cov.name] <- cov.ranges[2, cov.name]
      }
    }
  } else {
    X.ranges <- norm.cov.out$X.ranges
  }

  # For each covariate...
  for (cov.name in cont.cov.names) {

    # Get covariate of this iteration, alongside its min/max value
    if (!is.null(data)) {
      X <- data[, cov.name]
    }
    min.X <- X.ranges[1, cov.name]
    max.X <- X.ranges[2, cov.name]

    # Normalize the covariate based on the functions defined in the supplement
    # of Andrews, Shi (2013).
    if ((min.X > -Inf) & (max.X < Inf)) {
      if (is.null(norm.cov.out)) {
        X.normalized <- (X - min.X)/(max.X - min.X)
      }
      if (!is.null(x)) {
        x.elem <- x[cov.name]
        x.c.norm <- c(x.c.norm, (x.elem - min.X)/(max.X - min.X))
      }
    } else if ((min.X > -Inf) & (max.X == Inf)) {
      if (is.null(norm.cov.out)) {
        X.normalized <- (exp(X - min.X) - 1)/(1 + exp(X - min.X))
      }
      if (!is.null(x)) {
        x.elem <- x[cov.name]
        x.c.norm <- c(x.c.norm, (exp(x.elem - min.X) - 1)/(1 + exp(x.elem - min.X)))
      }
    } else if ((min.X == -Inf) & max.X < Inf) {
      if (is.null(norm.cov.out)) {
        X.normalized <- (2*exp(X - max.X))/(1 + exp(X - max.X))
      }
      if (!is.null(x)) {
        x.elem <- x[cov.name]
        x.c.norm <- c(x.c.norm, (2*exp(x.elem - max.X))/(1 + exp(x.elem - max.X)))
      }
    } else {
      if (is.null(norm.cov.out)) {
        X.normalized <- exp(X)/(1 - exp(X))
      }
      if (!is.null(x)) {
        x.elem <- x[cov.name]
        x.c.norm <- c(x.c.norm, exp(x.elem)/(1 - exp(x.elem)))
      }
    }

    # Transform x.norm to a named vector (in the first iteration, it will be a
    # list).
    x.c.norm <- unlist(x.c.norm)

    # Store the result
    if (is.null(norm.cov.out)) {
      normalized.data[, cov.name] <- X.normalized
    }
  }

  # If a previous output of this function was supplied, access the precomputed
  # normalized data
  if (!is.null(norm.cov.out)) {
    normalized.data <- norm.cov.out$normalized.data
  }

  # If x fell outside of the range of the covariates, the normalized values for
  # x might be negative. Note however that during the execution of the test by
  # Bei (2024), this should never occur.
  x.c.norm <- pmax(pmin(x.c.norm, 1), 0)

  # Reconstruct the full vector for x, with normalized continuous elements
  x.norm <- x
  x.norm[cont.cov.names] <- x.c.norm

  # Return the normalized data
  return(list("normalized.x" = x.norm,
              "normalized.data" = normalized.data,
              "X.ranges" = X.ranges,
              "covariate.names" = covariate.names,
              "cont.cov.names" = cont.cov.names))
}

#' @title Normalize the covariates of a data set to lie in the unit interval by
#' transforming based on PCA.
#'
#' @description This function normalized the covariates in the data to lie in
#' the unit interval based on a principal component analysis. It is useful to
#' perform this step when defining the instrumental functions later on. This
#' function is used in \code{G.box}, \code{G.spline} and by extension \code{G.cd}.
#'
#' @param data (optional) Data set to be used to construct the normalizing
#' transformation. Default is \code{data = NULL}.
#' @param x (optional) Vector of covariates to be normalized alongside the data.
#' Default is \code{x = NULL}.
#' @param idxs.c (optional) Vector of indices of covariates that are continuous.
#' Note that that indices are relative to the covariate vector, not the full
#' data set. Default value is \code{idxs.c = "all"}, which indicates that all
#' elements should be regarded as continuous. If \code{idxs.c = NULL}, all
#' elements are regarded as discrete.
#' @param norm.cov.out (optional) The output of a previous call to this function.
#' Can be used to speed up computation. If both \code{data} and
#' \code{norm.cov.out} are supplied to the function, the function will throw an
#' error. Default is \code{norm.cov.out = NULL}
#' @param ... Allows easier interchangeability between covariate normalization
#' functions. All arguments specified under \code{...} will be ignored.
#'
#' @import stats
#'
normalize.covariates2 <- function(data = NULL, x = NULL, idxs.c = "all",
                                  norm.cov.out = NULL, ...) {

  # Precondition checks
  if (is.null(data) & is.null(norm.cov.out)) {
    stop("Either data or norm.cov.out should be supplied to this function.")
  }
  if (!is.null(data) & !is.null(norm.cov.out)) {
    stop("Ambiguous function arguments: both data and norm.cov.out are supplied.")
  }

  # Extract the covariates from the data set, if applicable. Else extract the
  # necessary parameters from the previous function call.
  if (!is.null(data)) {

    # Get all covariates names
    cov.idxs <- which(grepl("X[[:digit:]]+", colnames(data)))[-1]
    covariate.names <- colnames(data)[cov.idxs]

    # Only retain names of continuous covariates
    idxs.c <- if (is(idxs.c, "character")) {1:length(covariate.names)} else {idxs.c}
    cont.cov.names <- covariate.names[idxs.c]

    # Obtain matrix of continuous covariates
    covariate.matrix <- as.matrix(data[, cont.cov.names, drop = FALSE])

  } else {
    covariate.names <- norm.cov.out$covariate.names
    cont.cov.names <- norm.cov.out$cont.cov.names
  }

  # If supplied, rename the entries of x.
  if (!is.null(x)) {
    names(x) <- covariate.names
  }

  # Define some useful variables
  n <- nrow(data)

  # Define function to transform unit circle towards [-1, 1]^d. Used to make
  # observations more evenly distributed in [-1, 1]^d later on. Else, extract
  # them from the previous output.
  transform.scores <- function(scores) {
    sin((1/2)*pi*scores)
  }

  # If the output of a previous call to this function was not supplied, obtain
  # the transformation parameters to be used.
  if (is.null(norm.cov.out)) {

    # If there are continuous covariates in the data set to be normalized...
    if (length(idxs.c > 0)) {

      # Apply principal component analysis to the covariates and obtain scores
      ev <- eigen(var(covariate.matrix))$vectors
      scores <- covariate.matrix %*% ev

      # Scale every covariate to have range of length 2
      scale <- apply(scores, 2, max) - apply(scores, 2, min)
      scale.mat <- matrix(rep(2/scale, n), ncol = ncol(scores), byrow = TRUE)
      scores.scaled <- scale.mat * scores

      # Shift the covariates into [-1, 1]^d
      shift <- apply(scores.scaled, 2, min) + 1
      shift.mat <- matrix(rep(shift, n), ncol = ncol(scores.scaled), byrow = TRUE)
      scores.scaled.shifted <- scores.scaled - shift.mat

      # Transform covariates to be more uniformly distributed in [-1, 1]^d
      scores.transformed <- transform.scores(scores.scaled.shifted)

      # Shift and scale the covariates into [0, 1]^d
      covariates.norm <- (1/2)*scores.transformed + (1/2)

      # Create the data frame with normalized covariates
      data.norm <- data
      data.norm[, cont.cov.names] <- covariates.norm

      # If there are no continuous covariates in the data to be normalized...
    } else {

      ev <- NULL
      scale <- NULL
      shift <- NULL
      data.norm <- data
    }

    # If the argument for norm.cov.out was provided...
  } else {

    ev <- norm.cov.out$ev
    scale <- norm.cov.out$scale
    shift <- norm.cov.out$shift
    data.norm <- norm.cov.out$normalized.data

  }

  # If x was supplied, transform x in the same way
  x.norm <- NULL
  if (!is.null(x)) {

    if (length(idxs.c) > 0) {
      # Transform the continuous elements in x
      x.c <- x[cont.cov.names]
      x.c.norm <- (1/2)*transform.scores(((matrix(x.c, nrow = 1) %*% ev) * (2/scale)) - shift) + (1/2)
      x.c.norm <- as.numeric(x.c.norm)

      # Construct entire vector
      x.norm <- x
      x.norm[cont.cov.names] <- x.c.norm
    } else {
      x.norm <- x
    }
  }

  # Return the results
  return(list("normalized.x" = x.norm,
              "normalized.data" = data.norm,
              "ev" = ev,
              "scale" = scale,
              "shift" = shift,
              "covariate.names" = covariate.names,
              "cont.cov.names" = cont.cov.names))
}

#' @title Get anchor points on which to base the instrumental functions
#'
#' @description The points returned by this function can be used as corner
#' points in the family of box functions, or as knots in the family of B-spline
#' functions.
#'
#' @param data Data set.
#' @param n.if.per.cov Number of instrumental functions to use per continuous
#' covariate.
#' @param normalized Boolean value indicating whether the covariates in the
#' given data frame have been normalized. Default is \code{normalized = FALSE}.
#'
get.anchor.points <- function(data, n.if.per.cov, normalized = FALSE) {

  # Get column indices of covariates in the data (excluding the intercept)
  cov.idxs <- which(grepl("X[[:digit:]]+", colnames(data)))[-1]

  # Get the number of covariates
  n.cov <- length(cov.idxs)

  # If the data should have normalized covariates, check that it is the case
  if (normalized) {
    for (cov.idx in cov.idxs) {
      if (!all((0 <= data[, cov.idx]) & (data[, cov.idx] <= 1))) {
        stop("Unnormalized data detected")
      }
    }
  }

  # Initialize object that will store the anchor points
  ap <- matrix(nrow = n.cov, ncol = n.if.per.cov + 1)
  rownames(ap) <- paste0("X", 1:n.cov)

  # For each covariate, determine an appropriate range for the boxes
  for (idx in 1:n.cov) {
    if (normalized) {

      # When the data is normalized, anchor points are evenly spaced in [0, 1]
      ap[idx, ] <- seq(0, 1, length.out = n.if.per.cov + 1)

    } else {

      # Select covariate of this iteration
      X <- data[, cov.idxs[idx]]

      # Determine the grid of corner points corresponding to this covariate for the
      # boxes.
      ap[idx, ] <- seq(min(X), max(X), length.out = n.if.per.cov + 1)

    }
  }

  # Return the results
  ap
}

#' @title Family of box functions
#'
#' @description This function defined the class of box functions as defined in
#' Willems et al. (2024+).
#'
#' @param x Vector of covariates to be normalized alongside the data. Default is
#' \code{x = NULL}.
#' @param g.idx Index of the instrumental function, in \{1, ..., n.inst.func\}.
#' @param data Data frame.
#' @param n.box.per.cov Number of box functions to consider per continuous
#' covariate.
#' @param norm.func Function to be used to normalize the covariates.
#' @param cov.ranges Matrix of ranges of the covariates. Used for normalizing
#' the data to the unit interval before applying the instrumental functions.
#' Default is \code{cov.ranges = NULL}.
#' @param norm.cov.out Output of a preliminary call to the supplied covariate
#' normalization function.
#' @param ... Additional arguments will be ignored. Useful for allowing
#' compatibility with the implementations of other instrument function families.
#' Specifically, it allows to ignore the \code{degree} argument used in
#' 'G.spline.R' and 'G.cd.R'.
#'
#' @importFrom EnvStats base
#'
G.box <- function(x, g.idx, data, n.box.per.cov, norm.func, cov.ranges = NULL,
                  norm.cov.out = NULL, ...) {

  # Normalize the covariates to lie in the unit interval
  if (is.null(norm.cov.out)) {
    out.norm <- norm.func(data = data, x = x, cov.ranges = cov.ranges,
                          norm.cov.out = NULL)
  } else {
    out.norm <- norm.func(data = NULL, x = x, cov.ranges = cov.ranges,
                          norm.cov.out = norm.cov.out)
  }
  x.norm <- out.norm[["normalized.x"]]

  # Get column indices of covariates in the data (excluding the intercept)
  cov.idxs <- which(grepl("X[[:digit:]]+", colnames(data)))[-1]

  # Get the number of covariates
  n.cov <- length(cov.idxs)

  # Construct matrix of corner points of the boxes
  cp <- matrix(rep(seq(0, 1, length.out = n.box.per.cov + 1), n.cov),
               nrow = n.cov, byrow = TRUE)
  rownames(cp) <- paste0("X", 1:n.cov)

  # For each dimension, determine the range of the covariate corresponding to
  # this box function
  range.idxs <- base(g.idx - 1, base = n.box.per.cov, num.digits = n.cov) + 1
  ranges <- NULL
  for (idx in 1:n.cov) {
    ranges <- rbind(ranges, c(cp[idx, range.idxs[idx]], cp[idx, range.idxs[idx] + 1]))
  }

  # Return the results
  as.numeric(all(ranges[,1] <= x.norm) & all(x.norm <= ranges[,2]))
}

#' @title Evaluate the specified B-spline, defined on the unit interval
#'
#' @description This function evaluates the specified B-spline defined on the
#' unit interval, when considering \code{n.if.per.cov} B-splines. Currently, the
#' implementation is based on the one in Andrews, Shi 2013 (supplementary
#' materials).
#'
#' @param x value inside the unit interval at which to evaluate the spline.
#' @param spline.index Index of the spline to evaluate.
#' @param n.if.per.cov Number of B-splines to consider over the unit interval.
#' @param degree Degree of the B-splines. Default is \code{degree = 3}.
#'
#' @importFrom splines2 bSpline
#'
#' @references Andrews, D.W.K. and Shi, X. (2013). Inference based on
#' confitional moment inequalities. Econometrica. 81(2):609-666.
#'
Bspline.unit.interval <- function(x, spline.index, n.if.per.cov, degree = 3) {

  # Precondition checks
  if (n.if.per.cov <= degree) {
    stop("n.if.per.cov must be larger than degree")
  }

  # Create vector of (boundary) knots to be used to construct the spline
  knots <- seq(0, 1, length.out = n.if.per.cov + 1 - degree)
  width <- diff(knots)[1]
  if (degree > 1) {
    knots <- c(-((degree - 1):1) * width, knots, 1 + 1:(degree - 1) * width)
  }
  if (degree != 0) {
    boundary.knots <- c(-degree * width, 1 + degree * width)
  } else {
    boundary.knots <- c(-width, 1)
    knots <- knots[-c(length(knots))]
  }

  # Obtain all spline function evaluations
  spline.evals <- bSpline(x, knots = knots, Boundary.knots = boundary.knots,
                          degree = degree)

  # Subset spline functions to the ones inside the unit interval
  spline.evals <- spline.evals[max(1, degree):(length(spline.evals) - degree)]

  # Return the requested spline function
  spline.evals[spline.index]
}

#' @title Family of spline instrumental functions
#'
#' @description This function normalizes the covariates to lie in the unit
#' interval and then evaluates each B-spline at each observation, multiplying
#' together the results per observation.
#'
#' @param x The vector of covariates at which to evaluate the B-splines
#' @param g.idx The index of the instrumental function. Note that g.idx ranges
#' between 1 and n.if.per.cov^n.cov, as an instrumental function is the product
#' of the appropriate B-spline evaluation for each element in the covariate
#' vector.
#' @param data Data frame containing the data.
#' @param n.if.per.cov Number of instrumental variables to be used per covariate.
#' @param norm.func Function to be used to normalize the covariates.
#' @param cov.ranges Matrix of ranges of the covariates. Used for normalizing
#' the covariates. If \code{cov.ranges = NULL}, the data will be normalized in a
#' data-dependent way. Default is \code{cov.ranges = NULL}.
#' @param norm.cov.out Output of a preliminary call to the given covariate
#' normalization function. Default is \code{norm.cov.out = NULL}.
#' @param degree Degree of B-splines to use. Default value is \code{degree = 3}.
#'
#' @importFrom EnvStats base
#'
G.spline <- function(x, g.idx, data, n.if.per.cov, norm.func, cov.ranges = NULL,
                     norm.cov.out = NULL, degree = 3) {

  # Get the number of covariates
  n.cov <- sum(grepl("X[1-9][[:digit:]]*", colnames(data)))

  # Normalize the covariates to lie in the unit interval
  if (is.null(norm.cov.out)) {
    out.norm <- norm.func(data = data, x = x, cov.ranges = cov.ranges,
                          norm.cov.out = NULL)
  } else {
    out.norm <- norm.func(data = NULL, x = x, cov.ranges = cov.ranges,
                          norm.cov.out = norm.cov.out)
  }
  x.norm <- out.norm[["normalized.x"]]

  # Get the index of the B-spline to be used for each of the covariates
  spline.idxs <- base(g.idx - 1, base = n.if.per.cov, num.digits = n.cov) + 1

  # Evaluate each element of the normalized covariate vector on the appropriate
  # B-spline. Multiply all results.
  spline.args <- cbind(x.norm, spline.idxs, n.if.per.cov, degree)
  spline.wrapper <- function(inp) {Bspline.unit.interval(inp[1], inp[2], inp[3], inp[4])}
  spline.evals <- apply(X = spline.args, MARGIN = 1, FUN = spline.wrapper)

  # Return the results
  prod(spline.evals)
}

#' @title Family of continuous/discrete instrumental function
#'
#' @description The function normalizes the continuous covariates to lie in the
#' unit interval and then evaluates the subvector of continuous covariates on
#' the specified family of instrumental function. For the discrete elements,
#' indicator functions are used for each level.
#'
#' @param x The vector of covariates at which to evaluate the B-splines
#' @param g.idx The index of the instrumental function.
#' @param data Data frame containing the data.
#' @param n.if.per.cov Number of instrumental functions per continuous covariate.
#' @param idxs.c Vector of indices of the continuous elements in the vector of
#' covariates.
#' @param G.c Family of instrumental functions to use for the subvector of
#' continuous covariates.
#' @param norm.func Function to be used to normalize the covariates.
#' @param discrete.covariate.levels Matrix containing as rows all possible
#' 'combined' levels of the discrete covariates. Default is
#' \code{discrete.covariate.levels = NULL}.
#' @param cov.ranges Matrix containing as its rows the lower and upper bounds
#' for each continuous covariate. Default is \code{cov.ranges = NULL}.
#' @param norm.cov.out Output of a preliminary call to a covariate normalization
#' function (defined above). This is used to speed up computations. Note that
#' this argument should only apply to continuous covariates!! Default is
#' \code{norm.cov.out = NULL}.
#' @param degree Degree of the spline functions to be used as instrumental
#' functions for the continuous covariates (if applicable). Default is
#' \code{degree = 3}.
#'
G.cd <- function(x, g.idx, data, n.if.per.cov, idxs.c, G.c, norm.func,
                 discrete.covariate.levels = NULL, cov.ranges = NULL,
                 norm.cov.out = NULL, degree = 3) {

  # Get the number of covariates
  n.cov <- sum(grepl("X[1-9][[:digit:]]*", colnames(data)))

  # Obtain subvectors of continuous and discrete covariates
  x.c <- x[idxs.c]
  x.d <- x[setdiff(1:length(x), idxs.c)]

  # Subset data to continuous and discrete covariates (also catching the cases
  # where there are no continuous/discrete variables).
  names.cov.c <- setdiff(paste("X", idxs.c, sep = ""), "X")
  names.cov.d <- setdiff(paste("X", setdiff(1:n.cov, idxs.c), sep = ""), "X")
  data.c <- data[, !(colnames(data) %in% names.cov.d)]
  data.d <- data[, !(colnames(data) %in% names.cov.c)]

  # Obtain the indices for the discrete instrumental functions
  n.inst.func.d <- max(nrow(unique(data[, names.cov.d, drop = FALSE])), 1)
  g.idx.d <- ((g.idx - 1) %% n.inst.func.d) + 1

  # Obtain the indices for the continuous instrumental functions
  n.inst.func.c <- n.if.per.cov^length(x.c)
  g.idx.c <- ((g.idx - 1) %/% n.inst.func.d) + 1

  # Subset the matrix of covariate ranges to the continuous covariates
  cov.ranges.c <- cov.ranges[, colnames(cov.ranges) %in% names.cov.c, drop = FALSE]

  # Obtain the instrumental function evaluation for the continuous elements
  eval.c <- 1
  if (length(x.c) > 0) {
    if (length(norm.cov.out$covariate.names) == length(norm.cov.out$cont.cov.names)) {
      norm.cov.out.c <- norm.cov.out
    } else {
      stop("Provided norm.cov.out argument should only apply to the continuous cov.")
    }
    eval.c <- G.c(x.c, g.idx.c, data.c, n.if.per.cov, norm.func, degree = degree,
                  cov.ranges = cov.ranges.c, norm.cov.out = norm.cov.out.c)
  }

  # Obtain the instrumental function evaluations for the discrete elements
  eval.d <- 1
  if (length(x.d) > 0) {
    eval.d <- as.numeric(all(x.d == discrete.covariate.levels[g.idx.d,]))
  }

  # Return the result
  eval.c * eval.d
}

#' @title Family of discrete/continuous instrumental functions, in the case of
#' many covariates.
#'
#' @description This function defines the family of discrete/continuous
#' instrumental functions in the case of many covariates. It does so by
#' considering a instrumental functions for each pair of entries in the given
#' covariate vector.
#'
#' @param x The vector of covariates at which to evaluate the B-splines
#' @param g.idx The index of the instrumental function.
#' @param data Data frame containing the data.
#' @param n.if.per.cov Number of instrumental functions per continuous covariate.
#' @param idxs.c Vector of indices of the continuous elements in the vector of
#' covariates.
#' @param G.c Family of instrumental functions to use for the subvector of
#' continuous covariates.
#' @param norm.func Function to be used to normalize the covariates.
#' @param info.manycov Data frame containing some information about the global
#' structure of the instrumental functions of this class. If
#' \code{info.manycov = NULL}, it will be computed during execution. Default is
#' \code{info.manycov = NULL}.
#' @param cov.ranges Matrix containing as its rows the lower and upper bounds
#' for each continuous covariate. Default is \code{cov.ranges = NULL}.
#' @param norm.cov.out Output of function that normalizes the covariates.
#' @param degree Degree of the spline functions to be used as instrumental
#' functions for the continuous covariates (if applicable). Default is
#' \code{degree = 3}.
#' @param ... Arguments specified here will be ignored. Used for compatibility
#' with other instrumental function classes.
#'
G.cd.mc <- function(x, g.idx, data, n.if.per.cov, idxs.c, G.c, norm.func,
                    info.manycov = NULL, cov.ranges = NULL,
                    norm.cov.out = NULL, degree = 3, ...) {

  #### Precompute/preset some necessary variables ####

  # Obtain vector of covariate names
  cov.names <- colnames(data)[grep("X[1-9][[:digit:]]*$", colnames(data))]

  # If the necessary information is not pre-supplied...
  if (is.null(info.manycov)) {

    # Obtain each pair of covariates in the data. For each, determine the amount
    # of instrumental functions.
    info.manycov <- data.frame(cov.pair = character(), n.if = numeric())
    for (cov.name.idx1 in 1:length(cov.names)) {
      for (cov.name.idx2 in 2:length(cov.names)) {
        if (cov.name.idx2 > cov.name.idx1) {

          # Name of covariates in the pair
          cov.name1 <- cov.names[cov.name.idx1]
          cov.name2 <- cov.names[cov.name.idx2]

          # Number of instrumental functions for each
          n.if1 <- ifelse(cov.name.idx1 %in% idxs.c, n.if.per.cov, length(unique(data[, cov.name1])))
          n.if2 <- ifelse(cov.name.idx2 %in% idxs.c, n.if.per.cov, length(unique(data[, cov.name2])))

          # Total number of instrumental functions
          n.if <- n.if1 * n.if2

          # Add to information data frame
          row <- list(cov.pair = sprintf("%s, %s", cov.name1, cov.name2),
                      n.if = n.if)
          info.manycov <- rbind(info.manycov, row)
        }
      }
    }

    # Add supplementary rows and columns
    info.manycov <- cbind(idx = 1:nrow(info.manycov),
                          info.manycov,
                          cumsum = cumsum(info.manycov$n.if))
    info.manycov <- rbind(list(idx = 0, cov.pair = "", n.if = 0, cumsum = 0),
                          info.manycov)
  }

  #### Select the relevant variables ####

  # Get pair of covariates corresponding to the given index of instrumental
  # function.
  cov.pair <- info.manycov[min(which(info.manycov$cumsum >= g.idx)), "cov.pair"]
  pair.vec <- strsplit(cov.pair, split = ", ")[[1]]

  # Get subset of data and covariate vector corresponding to the covariates
  data.sub <- data[, c("Y", "Delta", "X0", pair.vec)]
  x.sub <- x[which(cov.names %in% pair.vec)]
  g.idx.sub <- g.idx - info.manycov[min(which(info.manycov$cumsum >= g.idx)) - 1, "cumsum"]
  cov.ranges.sub <- cov.ranges[, pair.vec]

  #### Construct the class of instrumental functions for this pair ####

  # If both variables in the pair are continuous...
  if (all(which(cov.names %in% pair.vec) %in% idxs.c)) {
    eval <- G.c(x = x.sub, g.idx = g.idx.sub, data = data.sub,
                n.if.per.cov = n.if.per.cov, norm.func = norm.func,
                cov.ranges = cov.ranges.sub, degree = degree)

    # If both variables in the pair are binary...
  } else if (all(!(which(cov.names %in% pair.vec) %in% idxs.c))) {
    levels.var1 <- unique(data[, pair.vec[1]])
    levels.var2 <- unique(data[, pair.vec[2]])
    discrete.covariate.levels <- expand.grid(levels.var1, levels.var2)
    discrete.covariate.levels[order(discrete.covariate.levels[,1],
                                    discrete.covariate.levels[,2],
                                    decreasing = TRUE),]
    eval <- as.numeric(all(discrete.covariate.levels[g.idx.sub, ] == x.sub))

    # If one variable in the pair in continuous and the other one is binary...
  } else {
    cont.var.idx <- which(pair.vec %in% cov.names[idxs.c])
    disc.var.idx <- setdiff(1:2, cont.var.idx)
    levels.disc <- sort(unique(data[, pair.vec[disc.var.idx]]))
    n.levels.disc <- length(levels.disc)

    g.idx.sub.d <- ((g.idx.sub - 1) %% n.levels.disc) + 1
    g.idx.sub.c <- ((g.idx.sub - 1) %/% n.levels.disc) + 1
    data.sub.c <- data[, c("Y", "Delta", "X0", pair.vec[cont.var.idx])]
    cov.ranges.sub.c <- cov.ranges.sub[, pair.vec[cont.var.idx]]

    eval.d <- as.numeric(x.sub[disc.var.idx] == levels.disc[g.idx.sub.d])
    eval.c <- G.c(x = x.sub[cont.var.idx], g.idx = g.idx.sub.c,
                  data = data.sub.c, n.if.per.cov = n.if.per.cov,
                  norm.func = norm.func, cov.ranges = cov.ranges.sub.c,
                  degree = degree)

    eval <- eval.d * eval.c
  }

  #### Return the result ####

  eval
}

#' @title Evaluate each instrumental function at each of the observations.
#'
#' @description Obtain the evaluations of each observation on each of the
#' instrumental functions. (Used in function get.mi.mat.R)
#'
#' @param data Data frame.
#' @param hp List of hyperparameters. Notably, it contains the instrumental
#' function to be used in an element named \code{G}.
#'
get.instrumental.function.evals <- function(data, hp) {

  # Unpack hyperparameters
  n.inst.func <- hp[["n.inst.func"]]
  G <- hp[["G"]]

  # Initialize matrix that will store the instrumental function evaluations
  inst.func.evals <- matrix(nrow = nrow(data), ncol = n.inst.func)

  for (i in 1:nrow(data)) {

    # Get the covariate values of the i-th observation. Leave out the intercept.
    X <- as.matrix(data[i, grepl("X[[:digit:]]+", colnames(data))])
    X.no_int <- X[-1]

    # For each instrumental function, evaluate it at the covariates values of
    # the i-th observation.
    for (j in 1:n.inst.func) {
      inst.func.evals[i, j] <- G(X.no_int, j)
    }
  }

  # Return the results
  inst.func.evals
}

#### Low level functions: moment functions + derivatives ####

#' @title [DEPRECATED] Component function of the vector of moment functions m.
#'
#' @description THIS FUNCTION IS DEPRECATED AND WILL THROW A WARNING WHEN USED.
#' Use the faster function 'get.mi.mat.R' instead
#'
#' @param i Index of observation
#' @param j Index of moment function.
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the moment function
#' @param hp List of hyperparamerers.
#'
#' @noRd
#'
m.comp <- function(i, j, data, beta, t, hp) {

  # This function is deprecated
  warning("Attempted to use a deprecated function (m.comp.R)")

  # Unpack data
  Y <- data[i, "Y"]
  Delta <- data[i, "Delta"]
  X <- as.matrix(data[i, grepl("X[[:digit:]]+", colnames(data))])
  X.no_int <- matrix(X[-1], nrow = 1)

  # Unpack hyperparameters
  Lambda <- hp[["Lambda"]]
  G <- hp[["G"]]
  n.inst.func <- hp[["n.inst.func"]]

  # Compute moment function
  if (j <= n.inst.func) {
    (Lambda(X %*% beta) - as.numeric(Y <= t & Delta == 1)) * G(X.no_int, j)
  } else {
    (as.numeric(Y <= t) - Lambda(X %*% beta) ) * G(X.no_int, j - n.inst.func)
  }
}

#' @title [DEPRECATED] Vector of moment functions
#'
#' @description THIS FUNCTION IS DEPRECATED AND WILL THROW A WARNING WHEN USED.
#' Use the faster function 'get.mi.mat.R' instead
#'
#' @param i Index of observation
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the moment function.
#'
#' @noRd
#'
m <- function(i, data, beta, t, hp) {

  # This function is deprecated
  warning("Attempted to use a deprecated function (m.R).")

  # Number of instrumental functions
  n.inst.func <- hp[["n.inst.func"]]

  # Create vector of moment functions evaluated at (data, theta)
  rtrn <- c()
  for (j in 1:(2*n.inst.func)) {
    rtrn <- c(rtrn, m.comp(i, j, data, beta, t, hp))
  }

  # Return the results
  rtrn
}

#' @title Compute the conditional moment evaluations
#'
#' @description This function computes the 1(Y <= t) - Lambda(X^T beta(t)) and
#' Lambda(X^T beta(t)) - 1(Y <= t, Delta = 1) parts of the moment functions.
#' (Used in function get.mi.mat.R)
#'
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point of interest.
#' @param hp List of hyperparameters.
#'
#' @returns A vector of 2n elements containing in the first n positions the
#' evaluations of 1(Y <= t) - Lambda(X^T beta(t)) and in the last n positions
#' the evaluations of Lambda(X^T beta(t)) - 1(Y <= t, Delta = 1).
#'
get.cond.moment.evals <- function(data, beta, t, hp) {

  # Unpack hyperparameters
  Lambda <- hp[["Lambda"]]

  # Initialize matrix that will store the results
  evals <- matrix(nrow = nrow(data), ncol = 2)

  # For each observation, compute the evaluation of the two conditional moment
  # functions
  for (i in 1:nrow(data)) {

    # Get the values pertaining to the i-th observation
    Y <- data[i, "Y"]
    Delta <- data[i, "Delta"]
    X <- as.matrix(data[i, grepl("X[[:digit:]]+", colnames(data))])

    # Compute moment functions
    evals[i, 1] <- Lambda(X %*% beta) - as.numeric(Y <= t & Delta == 1)
    evals[i, 2] <- as.numeric(Y <= t) - Lambda(X %*% beta)
  }

  # Return the results
  evals
}

#' @title Faster implementation of vector of moment functions.
#'
#' @description
#' This function obtains the moment function evaluations.
#'
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the moment function. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters
#' @param inst.func.evals Matrix of instrumental function evaluations. If
#' \code{NULL}, it will be computed during execution of this function. Default
#' value is \code{inst.func.evals = NULL}.
#'
get.mi.mat <- function(data, beta, t, hp, inst.func.evals = NULL) {

  # Extract hyperparameters
  n.inst.func <- hp[["n.inst.func"]]
  n.cov <- sum(grepl("[1-9][[:digit:]]*$", colnames(data)))

  # Get instrumental function evaluations
  if (is.null(inst.func.evals)) {
    inst.func.evals <- t(get.instrumental.function.evals(data, hp))
  }

  # Create matrix of replicates of instrumental function evaluation. I.e.
  # ife = [inst.func.evals
  #        inst.func.evals
  #        ...
  #        inst.func.evals].
  ife <- do.call(rbind, replicate(2*length(t), inst.func.evals, simplify=FALSE))

  # Get conditional moment evaluations at each time point
  cmfe <- NULL
  for (time.point in t) {
    if (is(beta, "function")) {
      beta.t <- beta(time.point)
    } else if (length(t) == 1) {
      beta.t <- beta
    } else {
      beta.t <- beta[c(which(t == time.point), (length(t) + 1):length(beta))]
    }
    cond.m.evals <- get.cond.moment.evals(data, beta.t, time.point, hp)
    cmfe1 <- matrix(rep(cond.m.evals[,1], n.inst.func), nrow = n.inst.func, byrow = TRUE)
    cmfe2 <- matrix(rep(cond.m.evals[,2], n.inst.func), nrow = n.inst.func, byrow = TRUE)
    cmfe <- rbind(cmfe, cmfe1, cmfe2)
  }

  # Combine the conditional moment function with the instrumental function
  # evaluations and return the result.
  cmfe * ife
}

#' @title Vector of sample average of each moment function
#' \eqn{(\bar{m}_n(\theta))}.
#'
#' @description This function obtains the vector of sample averages of each
#' moment function.
#'
#'
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the moment functions. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param mi.mat Matrix of moment function evaluations. Can be used to avoid
#' some computation. Default is \code{mi.mat = NULL}.
#'
m.bar <- function(data, beta, t, hp, mi.mat = NULL) {

  # Number of instrumental functions
  n.inst.func <- hp[["n.inst.func"]]

  # Initialize vector that will contain the sum of all moment function
  # evaluations.
  m.evals.sum <- rep(0, 2*length(t)*n.inst.func)

  # Obtain the sum
  for (i in 1:nrow(data)) {
    if (is.null(mi.mat)) {
      mi <- NULL
      for (time.point in t) {
        mi <- c(mi, m(i, data, beta, time.point, hp))
      }
    } else {
      mi <- mi.mat[,i]
    }
    m.evals.sum <- m.evals.sum + mi
  }

  # Return the average
  m.evals.sum/nrow(data)
}

#' @title [DEPRECATED] Component function of the vector of derivatives of moment
#' functions (with respect to \eqn{\beta}).
#'
#' @description This function obtains the vector of partial derivatives of a
#' moment function, evaluated at a specified observation. This function is
#' deprecated.
#'
#'
#' @param i Index of observation
#' @param j Index of moment function.
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the derivative of the moment function
#' (not actually used in the implementation below)
#' @param hp List of hyperparameters
#'
#' @returns A vector containing the partial derivatives of the selected moment
#' function, evaluated at the specified observation.
#'
#' @noRd
#'
dm.comp <- function(i, j, data, beta, t, hp) {

  # Unpack data
  Y <- data[i, "Y"]
  Delta <- data[i, "Delta"]
  X <- as.matrix(data[i, grepl("X[[:digit:]]+", colnames(data))])
  X.no_int <- matrix(X[-1], nrow = 1)

  # Unpack hyperparameters
  dLambda <- hp[["dLambda"]]
  G <- hp[["G"]]
  n.inst.func <- hp[["n.inst.func"]]

  # Compute vector derivative of moment function
  if (j <= n.inst.func) {
    G(X.no_int, j) * dLambda(as.numeric(X %*% beta)) * X
  } else {
    - G(X.no_int, j - n.inst.func) * dLambda(as.numeric(X %*% beta)) * X
  }
}

#' @title [DEPRECATED] Vector of derivatives of moment functions
#'
#' @description This function returns a matrix containing the partial
#' derivatives of each moment function, evaluated at the specified observation.
#'
#' @param i Index of observation
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the derivative of the moment function
#' @param hp List of hyperparameters.
#'
#' @returns A matrix containing the partial derivatives of each moment
#' function, evaluated at the specified observation. Each row corresponds to a
#' moment function, each column corresponds to a coefficient.
#'
#' @noRd
#'
dm <- function(i, data, beta, t, hp) {

  # Warn user that the function is deprecated and might not work
  warning("Using deprecated function 'dm'!")

  # Number of instrumental functions
  n.inst.func <- hp[["n.inst.func"]]

  # Create vector of moment functions evaluated at (data, theta)
  rtrn <- NULL
  for (j in 1:(2*n.inst.func)) {
    rtrn <- rbind(rtrn, dm.comp(i, j, data, beta, t, hp))
  }

  # Return the results
  rtrn
}

#' @title Matrix of derivatives of conditional moment functions
#'
#' @description This function evaluates the derivatives of the conditional
#' moment function at each observation. Used in get.dmi.tens.R
#'
#' @param data Data frame.
#' @param beta Parameter vector.
#' @param t Time point of interest.
#' @param hp List of hyperparameters.
#'
get.deriv.mom.func <- function(data, beta, t, hp) {

  # Extract hyperparameters
  n.inst.func <- hp[["n.inst.func"]]
  dLambda <- hp[["dLambda"]]
  n.param <- length(beta)
  n <- nrow(data)

  # Initialize matrix that will store the results
  evals.m1 <- matrix(nrow = n, ncol = n.param)
  evals.m2 <- matrix(nrow = n, ncol = n.param)

  # For each observation, compute the evaluation of the two conditional moment
  # functions
  for (i in 1:nrow(data)) {

    # Get the values pertaining to the i-th observation
    Y <- data[i, "Y"]
    Delta <- data[i, "Delta"]
    X <- as.matrix(data[i, grepl("X[[:digit:]]+", colnames(data))])

    # Compute derivatives of moment functions
    evals.m1[i,] <- dLambda(as.numeric(X %*% beta)) * X
  }
  evals.m2 <- -evals.m1

  # Return the results
  evals <- array(dim = c(n, n.param, 2))
  evals[, , 1] <- evals.m1
  evals[, , 2] <- evals.m2
  evals
}

#' @title Faster implementation to obtain the tensor of the evaluations of the
#' derivatives of the moment functions at each observation.
#'
#' @description This function provides a faster implementation of obtaining the
#' evaluations of the derivative of the moment functions at each observation
#' (wrt the previous implementation using 'dm.comp' and 'dm.R'). Used in the
#' function G.hat.R
#'
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point of interest. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param inst.func.evals Precomputed matrix of instrumental function
#' evaluations. Defaults is \code{inst.func.evals = NULL}, in which case the
#' evaluations will be done inside this function.
#'
get.dmi.tens <- function(data, beta, t, hp, inst.func.evals = NULL) {

  # Extract hyperparameters
  n.inst.func <- hp[["n.inst.func"]]
  n.param <- length(beta)
  n <- nrow(data)

  # Compute the instrumental function evaluations if necessary
  if (is.null(inst.func.evals)) {
    inst.func.evals <- t(get.instrumental.function.evals(data, hp))
  }

  # Repeat the matrix of instrumental function evaluations into a tensor of the
  # correct dimension (note that the instrumental function evaluations do not
  # depend on the covariates).
  inst.func.tens <- array(dim = c(2*n.inst.func*length(t), n.param, n))
  for (i in 1:n.param) {
    inst.func.tens.i <- inst.func.evals
    for (dummy in 2:(2*length(t))) {
      inst.func.tens.i <- rbind(inst.func.tens.i, inst.func.evals)
    }
    inst.func.tens[, i, ] <- inst.func.tens.i
  }

  # Compute matrix of derivatives of conditional moment functions, evaluated at
  # each observation.
  deriv.mom.evals.list <- list()
  for (time.point in t) {
    if (is(beta, "function")) {
      beta.t <- beta(time.point)
    } else if (length(t) == 1) {
      beta.t <- beta
    } else {
      beta.t <- beta[c(which(t == time.point), (length(t)+1):length(beta))]
    }

    # Compute the derivatives of the moment functions at each time point wrt the
    # appropriate intercept parameter. Derivatives wrt intercept parameters at
    # other time points are zero
    deriv.mom.func <- array(0, dim = c(n, n.param, 2))
    deriv.mom.func[, c(which(t == time.point), (length(t)+1):n.param), ] <-
      get.deriv.mom.func(data, beta.t, time.point, hp)

    # Store the result
    deriv.mom.evals.list[[as.character(time.point)]] <- deriv.mom.func
  }

  # Create tensor of evaluations of the derivatives of the moment functions
  deriv.cond.mom.func <- array(dim = c(2*length(t)*n.inst.func, n.param, n))
  for (time.point.idx in 1:length(deriv.mom.evals.list)) {
    for (j in 1:n.inst.func) {
      deriv.cond.mom.func[2*n.inst.func*(time.point.idx - 1) + j, ,] <-
        t(deriv.mom.evals.list[[time.point.idx]][, , 1])
      deriv.cond.mom.func[2*n.inst.func*(time.point.idx - 1) + j + n.inst.func, ,] <-
        t(deriv.mom.evals.list[[time.point.idx]][, , 2])
    }
  }

  # Compute tensor of evaluations of derivatives of unconditional moment
  # functions.
  dmi.tens <- inst.func.tens * deriv.cond.mom.func
}

#' @title Vector of sample average of each moment function
#' \eqn{(\bar{m}_n(\theta))}.
#'
#' @description This function computes the matrix containing the sample average
#' of the partial derivatives of the moment functions.
#'
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the derivative of the moment
#' functions. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param dmi.tens Tensor of derivative moment function evaluations. Can be used
#' to avoid some computation. Default is \code{dmi.tens = NULL}.
#'
#' @returns A matrix containing the sample average of the partial derivatives of
#' the moment functions. Each row corresponds to a moment function, each column
#' corresponds to a coefficient.
dm.bar <- function(data, beta, t, hp, dmi.tens = NULL) {

  # Number of instrumental functions
  n.inst.func <- hp[["n.inst.func"]]

  # Number of covariates
  n.param <- length(beta)

  # Initialize vector that will contain the sum of all moment function
  # evaluations.
  dm.evals.sum <- matrix(0, nrow = 2*length(t)*n.inst.func, ncol = n.param)

  # Obtain the sum
  for (i in 1:nrow(data)) {
    if (is.null(dmi.tens)) {
      dmi <- dm(i, data, beta, t, hp)
    } else {
      dmi <- dmi.tens[,,i]
    }
    dm.evals.sum <- dm.evals.sum + dmi
  }

  # Return the average
  dm.evals.sum/nrow(data)
}

#### Low level functions: variances/correlation of moment functions ####

#' @title Compute the variance-covariance matrix of the moment functions.
#'
#' @description This function comptutes the empricical variance-covariance
#' matrix of the moment functions.
#'
#' @param data Data frame.
#' @param beta Coefficient vector.
#' @param t Time point of interest.
#' @param hp List of hyperparameters.
#' @param m.avg A precomputed vector of the sample average of the moment
#' functions. If not supplied, this vector is computed. Default is
#' \code{m.avg = NULL}.
#' @param mi.mat A precomputed matrix of moment function evaluations at each
#' observation. If supplied, some computations can be skipped. Default is
#' \code{mi.mat = NULL}.
#'
Sigma.hat <- function(data, beta, t, hp, m.avg = NULL, mi.mat = NULL) {

  # Number of instrumental functions
  n.inst.func <- hp[["n.inst.func"]]

  # Sample average of the moment functions
  if (is.null(m.avg)) {
    m.avg <- m.bar(data, beta, t, hp)
  }

  # Initialize matrix that will contain the sum of all outer products used in
  # obtaining the sample variance-covariance matrix
  sig.evals.sum <- matrix(0, nrow = 2*length(t)*n.inst.func, ncol = 2*length(t)*n.inst.func)

  # Obtain the sum
  for (i in 1:nrow(data)) {
    if (is.null(mi.mat)) {
      mi <- m(i, data, beta, t, hp)
    } else {
      mi <- mi.mat[,i]
    }
    sig.evals.sum <- sig.evals.sum + outer(mi - m.avg, mi - m.avg)
  }

  # Return the average
  sig.evals.sum/nrow(data)
}

#' @title Obtain the diagonal matrix of sample variances of moment functions
#'
#' @description This function computes the diagonal matrix of the sample
#' variance-covariance matrix.
#'
#' @param input Can either be the variance-covariance matrix obtained from the
#' function Sigma.hat, or the data frame.
#' @param beta The coefficient vector. Only needs to be supplied when the
#' argument for \code{input} is the data frame.
#' @param t The time point of interest. Only needs to be supplied when the
#' argument for \code{input} is the data frame.
#' @param hp List of hyperparameters. Only needs to be supplied when the
#' argument for \code{input} is the data frame.
#' @param m.avg See documentation of \code{Sigma.hat}. Only needs to be supplied
#' when the argument for \code{input} is the data frama.
#' @param mi.mat See documentation of \code{Sigma.hat}. Only needs to be supplied
#' when the argument for \code{input} is the data frama.
#'
D.hat <- function(input, beta = NULL, t = NULL, hp = NULL, m.avg = NULL,
                  mi.mat = NULL) {

  # If the given input is a variance-covariance matrix...
  if (is.matrix(input)) {
    Sigma <- input
    return(diag(diag(Sigma), nrow = nrow(Sigma)))

    # If the given input is a data frame...
  } else {

    data <- input

    # beta, t and hp should be specified in this case.
    if (any(is.null(c(beta, t, hp)))) {
      stop("When input is not a matrix, beta, t and hp should be specified.")
    }

    # Evaluations of the moment functions at each observation
    if (is.null(mi.mat)) {
      for (i in 1:nrow(data)) {
        mi.mat <- cbind(mi.mat, m(i, data, beta, t, hp))
      }
    }

    # Sample average of the moment functions
    if (is.null(m.avg)) {
      m.avg <- m.bar(data, beta, t, hp, mi.mat = mi.mat)
    }

    # Initialize variable that will store the result
    sum.sq.dev <- rep(0, length(m.avg))

    # Obtain the sum of squared deviations from the sample mean
    for (i in 1:nrow(data)) {
      mi <- mi.mat[,i]
      sum.sq.dev <- sum.sq.dev + (mi - m.avg)^2
    }

    # Return the resuls
    return(diag(sum.sq.dev/nrow(data), nrow = length(m.avg)))
  }
}

#' @title Obtain the correlation matrix of the moment functions
#'
#' @description This function computes the correlation matrix corresponding to
#' the variance-covariance matrix as returned by \code{Sigma.hat.R}
#'
#' @param Sigma The output of the function Sigma.hat
#'
Omega.hat <- function(Sigma) {
  sqrt.D.hat.inv <- solve(sqrt(D.hat(Sigma)))
  sqrt.D.hat.inv %*% Sigma %*% sqrt.D.hat.inv
}

#' @title Obtain the matrix of partial derivatives of the sample variances.
#'
#' @description This function computes the matrix of sample derivatives of the
#' sample variances.
#'
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to evaluate the (derivatives of) the moment
#' functions. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparamerers.
#' @param mi.mat A precomputed matrix of moment function evaluations at each
#' observation. If supplied, some computations can be skipped. Default is
#' \code{mi.mat = NULL}.
#' @param m.avg A precomputed vector of the sample average of the moment
#' functions. If not supplied, this vector is computed. Default is
#' \code{m.avg = NULL}.
#' @param dm.avg Matrix of precomputed sample averages of the derivatives of the
#' moment functions. Default is \code{dm.avg = NULL}.
#' @param dmi.tens 3D tensor of precomputed evaluations of the derivatives of
#' the moment functions. Default is \code{dmi.tens = NULL}.
#'
#' @returns A matrix containing the partial derivatives of the variances of the
#' moment functions. Each row corresponds to a moment function, each column
#' corresponds to a covariate.
dD.hat <- function(data, beta, t, hp, mi.mat = NULL, m.avg = NULL,
                   dm.avg = NULL, dmi.tens = NULL) {

  # Define some useful variables
  n <- nrow(data)
  n.param <- length(beta)
  n.inst.func <- hp[["n.inst.func"]]

  # Evaluations of the moment functions at each observation
  if (is.null(mi.mat)) {
    for (i in 1:nrow(data)) {
      mi.mat <- cbind(mi.mat, m(i, data, beta, t, hp))
    }
  }

  # Sample average of the moment functions
  if (is.null(m.avg)) {
    m.avg <- m.bar(data, beta, t, hp, mi.mat = mi.mat)
  }

  # Evaluations of the derivatives of moment functions at each observation
  if (is.null(dmi.tens)) {
    dmi.tens <- array(dim = c(2*n.inst.func, n.param, n))
    for (i in 1:n) {
      dmi.tens[, ,i] <- dm(i, data, beta, t, hp)
    }
  }

  # Sample average of the derivatives of moment functions
  if (is.null(dm.avg)) {
    dm.avg <- dm.bar(data, beta, t, hp, dmi.tens = dmi.tens)
  }

  # Compute dD.hat
  sum <- 0
  for (i in 1:n) {
    mi <- mi.mat[,i]
    dmi <- dmi.tens[,,i]
    sum <- sum + 2 * matrix(rep(mi - m.avg, n.param), ncol = n.param) * (dmi - dm.avg)
  }

  # Return the result
  sum/n
}

#' @title Compute the Gn matrix in step 3b of Bei (2024).
#'
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to evaluate the (derivatives of) the moment
#' functions.
#' @param hp List of hyperparamerers.
#' @param mi.mat A precomputed matrix of moment function evaluations at each
#' observation. If supplied, some computations can be skipped. Default is
#' \code{mi.mat = NULL}.
#' @param m.avg A precomputed vector of the sample average of the moment
#' functions. If not supplied, this vector is computed. Default is
#' \code{m.avg = NULL}.
#' @param dm.avg Matrix of precomputed sample averages of the derivatives of the
#' moment functions. Default is \code{dm.avg = NULL}.
#' @param dmi.tens 3D tensor of precomputed evaluations of the derivatives of
#' the moment functions. Default is \code{dmi.tens = NULL}.
#' @param D Diagonal of D-matrix.
#'
#' @returns A matrix containing the partial derivatives of the variances of the
#' moment functions. Each row corresponds to a moment function, each column
#' corresponds to a covariate.
#'
#' @references Bei, X. (2024). Local linearieation based subvector inference in
#' moment inequality models. Journal of Econometrics. 238:105549-
#'
G.hat <- function(data, beta, t, hp, mi.mat = NULL, m.avg = NULL,
                  dm.avg = NULL, dmi.tens = NULL, D = NULL) {

  # Define some useful variables
  n <- nrow(data)
  n.param <- length(beta)
  n.inst.func <- hp[["n.inst.func"]]

  # Evaluations of the moment functions at each observation
  if (is.null(mi.mat)) {
    mi.mat <- get.mi.mat(data, beta, t, hp)
  }

  # Sample average of the moment functions
  if (is.null(m.avg)) {
    m.avg <- m.bar(data, beta, t, hp, mi.mat = mi.mat)
  }

  # Evaluations of the derivatives of moment functions at each observation
  if (is.null(dmi.tens)) {
    dmi.tens <- get.dmi.tens(data, beta, t, hp)
  }

  # Sample average of the derivatives of moment functions
  if (is.null(dm.avg)) {
    dm.avg <- dm.bar(data, beta, t, hp, dmi.tens = dmi.tens)
  }

  # Compute the diagonal of Dn
  if (is.null(D)) {
    Dn <- diag(D.hat(data, beta, t, hp, m.avg = m.avg, mi.mat = mi.mat))
  } else {
    Dn <- diag(D)
  }

  # Compute the derivative of Dn
  dDn <- dD.hat(data, beta, t, hp, mi.mat = mi.mat, m.avg = m.avg,
                dm.avg = dm.avg, dmi.tens = dmi.tens)

  # Compute Gn
  a <- dm.avg * matrix(rep(sqrt(Dn), n.param), ncol = n.param)
  b <- matrix(rep(m.avg, n.param), ncol = n.param) * (1/2) * matrix(rep(Dn^(-1/2), n.param), ncol = n.param) * dDn

  (a - b)/matrix(rep(Dn, n.param), ncol = n.param)
}

#### Low level functions: S functions ####

#' @title S-function
#'
#' @description This function computes the loss function at a given point.
#'
#' @param m Vector of averages of moment functions.
#' @param Sigma Sample variance-covariance matrix of moment functions.
#'
#' @returns S(m, Sigma).
#'
S.func <- function(m, Sigma) {

  # Number of moment functions
  p <- length(m)

  # Initialize variable
  S <- 0
  for (j in 1:p) {
    S <- S + max(-m[j]/sqrt(Sigma[j, j]), 0)^2
  }

  # Return results
  S
}

#### Low level functions: estimating the test statistic ####

#' @title 'Loss function' of the test statistic.
#'
#' @description This function implements the loss function used in computing
#' the test statistic.
#'
#' @param beta.sub Subvector of coefficient vector.
#' @param data Data frame.
#' @param t Time point of interest. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param c Unit vector containing unity at the location of the parameter of
#' interest.
#' @param r Value of the parameter of interest that is tested.
#' @param inst.func.evals Pre-computed matrix of insturmental function
#' evaluations. If not supplied, it will be computed during execution of this
#' function.
#'
#' @returns S-functions evaluation for the specified parameter vector.
#'
lf.ts <- function(beta.sub, data, t, hp, c, r, inst.func.evals = NULL) {

  # Sample size
  n <- nrow(data)

  # Minimum variance (used for computational reasons)
  min.var <- hp[["min.var"]]

  # Make the completed parameter vector
  if (length(t) == 1) {
    beta <- rep(r, length(beta.sub) + 1)
    beta[which(c == 0)] <- beta.sub
  } else {
    beta <- function(time.point) {
      beta <- rep(r, length(beta.sub) + 1)
      beta[which(c == 0)] <- beta.sub
      beta[1:length(t)] <- cumsum(beta[1:length(t)])
      beta <- beta[-which(t != time.point)]
      beta
    }
  }

  # Matrix of moment function evaluations
  mi.mat <- get.mi.mat(data, beta, t, hp, inst.func.evals)

  # Sample average of the moment functions
  m.avg <- m.bar(data, beta, t, hp, mi.mat = mi.mat)

  # Sample variance-covariance matrix
  svc <- Sigma.hat(data, beta, t, hp, mi.mat = mi.mat, m.avg = m.avg)

  # Ensure the invertibility of the sample variance-covariance matrix
  svc <- svc + min.var * diag(ncol(svc))

  # Sample variance diagonal matrix
  D <- D.hat(svc)

  # S-function
  S.func(sqrt(n) * diag(D)^(-1/2) * m.avg, Omega.hat(svc))
}

#' @title Obtain the test statistic by minimizing the S-function over the
#' feasible region \eqn{\beta(r)}.
#'
#' @param beta.init Starting value of minimization algorithm.
#' @param data Data frame.
#' @param par.space Matrix containing the bounds on the parameter space.
#' @param t Time point at which to evaluate beta. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param c Projection vector
#' @param r hypothesised value of the projection.
#' @param inst.func.evals Matrix of precomputed instrumental function
#' evaluations for each observation in the data set. If \code{NULL}, the
#' evaluations will be computed during execution of this function. Default is
#' \code{inst.func.evals = NULL}.
#'
#' @returns A list containing the value of the test statistic and the parameter
#' at which this value was attained.
#'
#' @import stats
#' @importFrom nloptr nloptr
#'
get.test.statistic <- function(beta.init, data, par.space, t, hp, c, r,
                               inst.func.evals = NULL) {

  # Define some useful parameters
  n.param <- length(c)

  # If data.init represents the full vector, check whether it satisfies the
  # constraint and transform it into the unconstrained subvector
  if (length(beta.init) == n.param) {
    if (beta.init[which(c == 1)] != r) {
      stop("Given beta vector does not satisfy the constraint")
    }
    beta.init <- beta.init[which(c == 0)]
  }

  # Precompute the instrumental function evaluations
  if (is.null(inst.func.evals)) {
    inst.func.evals <- t(get.instrumental.function.evals(data, hp))
  }

  # Estimate the test statistic and the minimizer
  use.optim <- FALSE # Use stats::optim for optimization. Testing shows that
  # this may miss the global optimum, leading to over-
  # rejection.
  if (use.optim) {
    out <- optim(beta.init, lf.ts, data = data, t = t, hp = hp, c = c, r = r,
                 inst.func.evals = inst.func.evals,
                 method = "L-BFGS-B", lower = par.space[which(c == 0), 1],
                 upper = par.space[which(c == 0), 2], control = list(maxit = 200))

    Jnrh <- out$value
    beta.hat <- rep(r, length(c))
    beta.hat[which(c == 0)] <- out$par
  } else {
    out <- nloptr(beta.init,
                  eval_f = lf.ts,
                  lb = par.space[which(c == 0), 1],
                  ub = par.space[which(c == 0), 2],
                  opts = list("algorithm" = "NLOPT_LN_NEWUOA_BOUND",
                              "xtol_rel" = 1e-4,
                              "maxeval" = 1000),
                  data = data, t = t, hp = hp, c = c, r = r,
                  inst.func.evals = inst.func.evals)
    Jnrh <- out$objective
    beta.hat <- rep(r, length(c))
    beta.hat[which(c == 0)] <- out$solution
  }

  # Return the results
  return(list(Jnrh, beta.hat))
}

#### Low level functions: calculating the critical value of the test statistic ####

#' @title Loss function to compute Delta(beta).
#'
#' @description This function defines the loss function used in computing the
#' penalized local linear approximation of the test statistic in order to
#' construct the bootstrap distribution of the test statistic.
#'
#' @param Delta.sub Subvector of Delta.
#' @param vnb Bootstrapped stochastic process.
#' @param phi Moment selection functions.
#' @param Gn First-order approximation matrix.
#' @param Omegan Correlation matrix of sample moment functions.
#' @param beta Coefficient vector.
#' @param c Projection vector.
#' @param r Value of projected coefficient vector.
#' @param data Data frame.
#' @param par.space Matrix containing the bounds on the parameter space.
#' @param epsilon.n Parameter used in constructing the feasible region as in
#' Example 4.1 in Bei (2024). Not used in this function.
#' @param lambda.n Weight of penalty term.
#'
#' @returns Loss function evaluation evaluated at the given Delta.
#'
#' @references Bei, X. (2024). Local linearieation based subvector inference in
#' moment inequality models. Journal of Econometrics. 238:105549-
#'
lf.delta.beta1 <- function(Delta.sub, vnb, phi, Gn, Omegan, beta, c, r, data,
                           par.space, epsilon.n, lambda.n) {

  # Extract some useful parameters
  n <- nrow(data)

  # Make the completed parameter vector
  Delta <- rep(0, length(Delta.sub) + 1)
  Delta[which(c == 0)] <- Delta.sub

  # Value of the loss function
  S.func(vnb + phi + Gn %*% Delta, Omegan) + lambda.n/n * sum(Delta^2)

}

#' @title Compute the critical value of the test statistic.
#'
#' @description This function computes the critical value following the
#' algorithm of Section 4.3 in Bei (2024).
#'
#' @param BetaI.r Matrix containing in its columns the minimizers of the
#' S-function leading to the test statistic.
#' @param data Data frame.
#' @param t Time point of interest. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param c Projection vector.
#' @param r Result of projection of parameter vector onto \code{c}.
#' @param par.space Bounds on the parameter space.
#' @param inst.func.evals Matrix of precomputed instrumental function
#' evaluations for each observation in the data set. If \code{NULL}, the
#' evaluations will be computed during execution of this function. Default is
#' \code{inst.func.evals = NULL}.
#' @param alpha Confidence level.
#'
#' @returns The critical value for the test statistic.
#'
#' @import stats
#' @importFrom EnvStats base
#'
#' @references Bei, X. (2024). Local linearieation based subvector inference in
#' moment inequality models. Journal of Econometrics. 238:105549-
#'
get.cvLLn <- function(BetaI.r, data, t, hp, c, r, par.space,
                      inst.func.evals = NULL, alpha = 0.95) {

  # Define variables that will be useful throughout
  n <- nrow(data)
  J <- hp[["n.inst.func"]]*2
  n.beta <- ncol(BetaI.r)
  n.param <- nrow(BetaI.r)
  B <- hp[["B"]]
  kappa.n <- hp[["kappa.n"]]
  epsilon.n <- hp[["epsilon.n"]]
  lambda.n <- hp[["lambda.n"]]
  min.var <- hp[["min.var"]]

  # Precompute instrumental function evaluations
  if (is.null(inst.func.evals)) {
    inst.func.evals <- t(get.instrumental.function.evals(data, hp))
  }

  # Precompute moment function evaluations for all parameters in BetaI.r
  mi.tens <- array(dim = c(J*length(t), n, n.beta))
  for (beta.idx in 1:n.beta) {
    beta <- BetaI.r[, beta.idx]
    mi.tens[, , beta.idx] <- get.mi.mat(data, beta, t, hp, inst.func.evals)
  }

  # Precompute sample averages of moment functions for all beta in BetaI.r
  m.avg.mat <- matrix(nrow = J*length(t), ncol = n.beta)
  for (beta.idx in 1:n.beta) {
    beta <- BetaI.r[, beta.idx]
    m.avg.mat[, beta.idx] <- m.bar(data, beta, t, hp,
                                   mi.mat = mi.tens[, , beta.idx])
  }

  # Precompute variance-covariance matrix of moment functions for all beta.
  # Ensure the invertibility of each.
  Sigma.tens <- array(dim = c(J*length(t), J*length(t), n.beta))
  for (beta.idx in 1:n.beta) {
    beta <- BetaI.r[, beta.idx]
    Sigma.tens[, , beta.idx] <- Sigma.hat(data, beta, t, hp,
                                          m.avg = m.avg.mat[, beta.idx],
                                          mi.mat = mi.tens[ , , beta.idx])

    Sigma.tens[, , beta.idx] <- Sigma.tens[, , beta.idx] + min.var * diag(ncol(Sigma.tens[, , beta.idx]))
  }

  # Precompute the variance diagonal matrices for all beta in BetaI.r (stored
  # as vectors)
  D.diag.mat <- matrix(nrow = J*length(t), ncol = n.beta)
  for (beta.idx in 1:n.beta) {
    D.diag.mat[, beta.idx] <- diag(Sigma.tens[, , beta.idx])
  }

  # Precompute the square root of the inverse of the diagonal variance matrices
  D.inv.sqrt.diag.mat <- D.diag.mat^(-1/2)

  # Precompute phi(xi(beta)) of each beta in BetaI.r
  phi.mat <- matrix(nrow = J*length(t), ncol = n.beta)
  for (beta.idx in 1:n.beta) {
    beta <- BetaI.r[, beta.idx]
    phi <- rep(0, J)
    for (j in 1:J) {
      phi[j] <- max(sqrt(n) * m.avg.mat[j, beta.idx] * D.inv.sqrt.diag.mat[j, beta.idx] / kappa.n, 0)
    }
    phi.mat[, beta.idx] <- phi
  }

  # Precompute Gn(theta)
  Gn.tens <- array(dim = c(J*length(t), n.param, n.beta))
  for (beta.idx in 1:n.beta) {
    beta <- BetaI.r[, beta.idx]
    Gn.tens[, , beta.idx] <- G.hat(data, beta, t, hp,
                                   mi.mat = mi.tens[ , , beta.idx],
                                   m.avg = m.avg.mat[ , beta.idx],
                                   D = diag(D.diag.mat[, beta.idx], nrow  = J*length(t)))
  }

  # Initialize object that will store all bootstrapped test statistics
  JnLLb.r_vct <- NULL

  # For each bootstrap iteration, compute the bootstrapped test statistic
  for (b in 1:B) {

    # Simulate i.i.d standard normal random variables.
    zeta <- rnorm(n)

    # Initialize matrix to store all S-function evaluations needed to determine
    # the bootstrapped test statistic.
    S.evals <- NULL

    # Loop over all beta vectors inside BetaI.r and compute the corresponding
    # S-function evaluation.
    for (col.idx in 1:n.beta) {

      # Select the beta corresponding to this iteration
      beta <- BetaI.r[, col.idx]

      # Matrix of moment function evaluations
      mi.mat <- mi.tens[, , col.idx]

      # Sample average of the moment functions
      m.avg <- m.avg.mat[, col.idx]

      # Sample variance-covariance matrix
      Sigman <- Sigma.tens[, , col.idx]

      # Sample correlation matrix
      Omegan <- Omega.hat(Sigman)

      # (Inverse square root) sample variance matrix
      D.inv.sqrt <- diag(D.inv.sqrt.diag.mat[, col.idx], J*length(t))

      # Compute vnb
      vnb <- rep(0, J*length(t))
      for (i in 1:n) {
        vnb <- vnb + sqrt(1/n) * D.inv.sqrt %*% (mi.mat[,i] - m.avg) * zeta[i]
      }

      # phi(xi(theta))
      phi <- phi.mat[, col.idx]

      # Compute \hat{G}_n
      Gn <- Gn.tens[, , col.idx]

      # Evaluate the loss function at the origin. Store the result in a matrix
      delta.search <- matrix(c(rep(0, n.param), S.func(vnb + phi, Omegan)), ncol = n.param + 1)
      colnames(delta.search) <- c(paste0("X", 0:(n.param - 1)), "val")

      # Find the minimum of the loss function, starting from each corner point of
      # the feasible region.
      for (comb.nbr in 1:2^(n.param - 1)) {

        # Define the parameter bounds for the optimization
        optim.lb <- sqrt(n) * (par.space[which(c == 0), 1] + epsilon.n - beta[which(c == 0)])
        optim.ub <- sqrt(n) * (par.space[which(c == 0), 2] - epsilon.n - beta[which(c == 0)])

        # Get the corner point corresponding to this iteration
        comb <- base(comb.nbr - 1, 2, num.digits = n.param - 1)
        corner.point <- (1 - comb) * optim.lb + comb * optim.ub

        # Perform the optimization
        out <- optim(corner.point, lf.delta.beta1, vnb = vnb, phi = phi, Gn = Gn,
                     Omegan = Omegan, beta = beta, c = c, r = r, data = data,
                     par.space = par.space, epsilon.n = epsilon.n,
                     lambda.n = lambda.n, method = "L-BFGS-B",
                     lower = optim.lb, upper = optim.ub)

        # Extract the results
        val <- out$value
        solution <- rep(0, length(c))
        solution[which(c == 0)] <- out$par
        delta.search <- rbind(delta.search, c(solution, val))
      }

      # Obtain the 'global' minimum
      Delta.beta <- delta.search[which.min(delta.search[, "val"]), 1:n.param]

      # S function evaluation
      S.evals <- c(S.evals, S.func(vnb + phi + Gn %*% Delta.beta, Omegan))
    }

    # Compute the bootstrapped LL test statistic
    JnLLb.r_vct <- c(JnLLb.r_vct, min(S.evals))
  }

  # Obtain the (1 - \alpha)-quantile of the bootstrap distribution
  cvLLn <- quantile(JnLLb.r_vct, probs = alpha)
  cvLLn
}

#### Low level functions: analyzing the results ####

#' @title Obtain identified set based on results of main estimation algorithm.
#'
#' @description
#' Takes the results of the main estimation algorithm as input and outputs the
#' 1 dimensional identified set.
#'
#' @param test.results Results of main algorithm.
#'
#' @noRd
#'
get.identified.set <- function(test.results) {

  # If no feasible points were found, return [-\infty, \infty]
  if (length(which(test.results[, 2] <= test.results[, 3])) == 0) {
    return(c(-Inf, Inf))
  }

  # This is precisely step 6 in the algorithm described in Bei, 2024
  lb <- min(test.results[which(test.results[, 2] <= test.results[, 3]), 1])
  ub <- max(test.results[which(test.results[, 2] <= test.results[, 3]), 1])

  # Return results
  c(lb, ub)
}



#### Search strategies: helper functions ####

#' @title Obtain next point for feasible point search.
#'
#' @description
#' Function to obtain the next point to evaluate in the search for a feasible
#' point. This function evaluates the entire parameter space of the component of
#' theta as evenly as possible. Used in the initialization step
#' (feasible_point_search.R)
#'
#' @param evaluations Matrix of evaluations performed so far.
#' @param lb.theta Lower bound on the parameter of interest.
#' @param ub.theta Upper bound on the parameter of interest.
#'
#' @returns Next point in the feasible point search.
#'
get.next.point <- function(evaluations, lb.theta, ub.theta) {

  # Determine whether a next evaluation point is needed. I.e., whether there
  # already is a feasible point in the given set of evaluated points.
  if (any(evaluations[, "t.stat"] <= evaluations[, "crit.val"])) {
    return(list(theta.next = NULL, idx.after = NULL, stop = TRUE))
  }

  # If not, get the next point to evaluate.

  # To start, append the bounds of the theta values to the vector of evaluated
  # theta points, if necessary. Also store flag of whether bounds were added
  lb.added <- FALSE
  ub.added <- FALSE
  eval.points.extended <- evaluations[, "theta"]
  if (!(lb.theta %in% evaluations[, "theta"])) {
    eval.points.extended <- c(lb.theta, eval.points.extended)
    lb.added <- TRUE
  }
  if (!(ub.theta %in% evaluations[, "theta"])) {
    eval.points.extended <- c(eval.points.extended, ub.theta)
    ub.added <- TRUE
  }

  # Obtain the lengths of in the intervals between two consecutive points
  int.len <- diff(eval.points.extended)

  # Obtain the indices of the unique elements
  idx.uniques <- which.unique(matrix(int.len, nrow = length(int.len)), tol = 1e-10)

  # If all elements in int.len are the same (and hence length(idx.uniques) = 1),
  # return the midpoint between lb.theta and the smallest evaluated point.
  if (length(idx.uniques) == 1) {
    return(list(theta.next = sum(eval.points.extended[1:2])/2,
                idx.after = ifelse(lb.added, 0, 1),
                stop = FALSE))
  }

  # If not, length(idx.uniques) will necessarily be 2. Return the point obtained
  # by adding the smallest of the two interval lengths to the end point of the
  # last smallest interval.
  return(list(theta.next = eval.points.extended[idx.uniques[1] + 1] + int.len[1],
              idx.after = idx.uniques[1] + ifelse(lb.added, 0, 1),
              stop = FALSE))
}

#' @title Insert row into a matrix at a given row index
#'
#' @description
#' Used in initalization step (feasible_point_search.R).
#'
#' @param evaluations Matrix of violation function evaluations.
#' @param row Row (evaluations) to be added to the evaluation matrix.
#' @param idx.after Index of the row of \code{evaluations} after which the given
#' row should be placed.
#'
#' @returns Evaluation matrix.
#'
insert.row <- function(evaluations, row, idx.after) {

  if (idx.after == 0) {
    evaluations <- rbind(row, evaluations)
  } else if (idx.after == nrow(evaluations)) {
    evaluations <- rbind(evaluations, row)
  } else {
    evaluations <- rbind(evaluations[1:idx.after,], row, evaluations[(idx.after + 1):nrow(evaluations), ])
  }

  evaluations
}

#' @title Expected improvement
#'
#' @description Used in the M-step (M_step.R). Note: predict(fit.krige, ...)
#' has weird beheviour when making predictions for a single value in terms of
#' standard error. We work around this issue in this implementation.
#'
#' @param theta Vector of coefficients.
#' @param test.fun Test function (cf. \code{EstimationAlgorithmBei.R}).
#' @param fit.krige Fitted Kriging model.
#' @param theta.hash Tentative optimal value for theta, i.e., the largest or
#' smallest feasible value for theta (if dir = 1 or dir = -1, respectively). A
#' 'feasible value' is one that satisfies all moment restrictions.
#' @param dir Search direction. \code{dir = 1} corresponds to looking for an
#' upper bound. \code{dir = -1} corresponds to looking for a lower bound.
#'
#' @returns The expected improvement.
#'
#' @importFrom stats pnorm
#'
EI <- function(theta, test.fun, fit.krige, theta.hash, dir) {

  # Select another point for which to make the prediction (see note in the
  # description section of the documentation above).
  other.point <- theta + 1

  # Predicted (test statistic - critical value) based on kriging model
  pred.kriging <- predict(fit.krige, matrix(c(theta, other.point), nrow = 2))
  violation.theta <- pred.kriging[1, "y"]

  # Standard deviation of prediction
  sL.theta <- sqrt(pred.kriging[1, "s"])

  # Expected improvement
  dir * (theta - theta.hash) * (1 - pnorm(violation.theta/sL.theta))
}

#' @title Draw initial set of starting values for optimizing the expected
#' improvement.
#'
#' @description  Used in the M-step (get.starting.values.R). ToDo: Adapt this
#' code so as to also perform sample space contractions as in the MatLab
#' implementation of Bei (2024).
#'
#' @param theta.hash Tentative optimal value for theta, i.e., the largest or
#' smallest feasible value for theta (if dir = 1 or dir = -1, respectively). A
#' 'feasible value' is one that satisfies all moment restrictions.
#' @param dir Search direction. \code{dir = 1} corresponds to looking for an
#' upper bound. \code{dir = -1} corresponds to looking for a lower bound.
#' @param hyperparams List of hyperparameters.
#' @param EI.fun Function used to compute the expected improvement. See also
#' \code{EI}.
#'
#' @returns Initial set of starting values.
#'
#' @import stats
#'
#' @references Bei, X. (2024). Local linearieation based subvector inference in
#' moment inequality models. Journal of Econometrics. 238:105549-
#'
draw.sv.init <- function(theta.hash, dir, hyperparams, EI.fun) {

  #### Extract the necessary sampling hyperparameters ####

  # Minimum and maximum distance of sampled points from current best theta
  min.dist <- hyperparams[["min.dist"]]
  max.dist <- hyperparams[["max.dist"]]

  # Bounds of the parameter space for theta
  theta.lb <- hyperparams[["theta.lb"]]
  theta.ub <- hyperparams[["theta.ub"]]

  # Total number of sampled points required in initial drawing process
  nbr.init.sample.points <- hyperparams[["nbr.init.sample.points"]]

  # Number of points sampled per iteration in the initial drawing process
  nbr.points.per.iter.init <- hyperparams[["nbr.points.per.iter.init"]]

  # Total number of uniformly drawn points in the initial set of starting
  # values
  nbr.init.unif <- hyperparams[["nbr.init.unif"]]

  #### Draw points from a box ####

  # Initialize matrix that will store all the drawn points and their expected
  # improvement.
  init.draws <- matrix(nrow = 0, ncol = 2)
  colnames(init.draws) <- c("init.val", "EI")

  # Iterate until a prespecified number of points have been drawn or the
  # minimum distance of sampled points from theta.hash exceeds the maximum
  # distance (at the end of each iteration, the maximum distance is halved).
  while ((nrow(init.draws) < nbr.init.sample.points) & (min.dist < max.dist)) {

    # Compute bounds in which to draw points [See ToDo].
    if (dir == 1) {
      lb.init.sample.space <- theta.hash
      ub.init.sample.space <- min(theta.hash + max.dist, theta.ub)
    } else if (dir == -1) {
      lb.init.sample.space <- max(theta.hash - max.dist, theta.lb)
      ub.init.sample.space <- theta.hash
    }

    # Sample points in [lb.init.sample.space, ub.init.sample.space]
    draws <- runif(nbr.points.per.iter.init, min = lb.init.sample.space,
                   max = ub.init.sample.space)

    # For each draw, compute the expected improvement. Store it.
    for (draw in draws) {
      init.draws <- rbind(init.draws, c(draw, EI.fun(draw)))
    }

    # Only keep the draws with positive expected improvement
    keep.idxs <- which(init.draws[, "EI"] > 1e-10)
    init.draws <- init.draws[keep.idxs, , drop = FALSE]

    # Reduce the maximum distance from theta.hash
    max.dist <- max.dist / 2
  }

  # Only keep the top 'nbr.init.sample.points' points with highest expected
  # improvement, taking into account that init.draws can still be empty.
  if (nrow(init.draws) > 0) {
    init.draws <- init.draws[order(init.draws[,2], decreasing = TRUE), , drop = FALSE]
    init.draws <- init.draws[1:min(nbr.init.sample.points, nrow(init.draws)), , drop = FALSE]
  }

  #### Draw points uniformly from the parameter space ####

  # Draw from uniform distribution
  unif.draws <- runif(nbr.init.unif, min = theta.lb, max = theta.ub)

  # Set the expected improvement of the uniformly drawn points equal to zero.
  # (This might not be true but we do not require EI to be positive for these
  # points, so there is no need to compute it)
  unif.draws.EI <- rep(0, nbr.init.unif)

  # Append to other starting values
  init.draws <- rbind(init.draws, cbind(unif.draws, unif.draws.EI))

  #### Return the results ####

  init.draws
}

#' @title Analogue to KMS_AUX4_MSpoints(...) in MATLAB code of Bei (2024).
#'
#' @description Create starting values for EI maximization. Used in the M-step
#' (get.starting.values.R).
#'
#' @param draws.init Initial draws.
#'
#' @references Bei, X. (2024). Local linearieation based subvector inference in
#' moment inequality models. Journal of Econometrics. 238:105549-
#'
MSpoint <- function(draws.init) {
  X <- draws.init[, 1]
  n.eval <- length(X)

  # Initialize object that will store ms_points
  ms_points <- matrix(0, nrow = n.eval, ncol = 2)

  # Sort the starting values
  temp <- sort(X)

  # Replace the k-th coordinate with averages with nearest neighbours
  temp1 <- c((temp[1] + temp[2])/2, (temp[-1] + temp[-length(temp)])/2)
  temp2 <- c(temp1[-1], (temp[length(temp) - 1] + temp[length(temp)])/2)
  sel <- rbinom(n.eval, 1, 0.5)

  # Randomize between candidates
  ms_points[ , 1] <- sel * temp1 + (1 - sel) * temp2

  # Use the ones that were not selected for the second dimension
  ms_points[ , 2] <- (1 - sel) * temp1 + sel * temp2

  # Reshape and return the evaluation points
  c(ms_points[, 1], ms_points[, 2])
}

#' @title Main function for obtaining the starting values of the expected
#' improvement maximization step.
#'
#' @description Obtain starting values used in the M-step (M_step.R).
#'
#' @param theta.hash Tentative optimal value for theta, i.e., the largest or
#' smallest feasible value for theta (if dir = 1 or dir = -1, respectively). A
#' 'feasible value' is one that satisfies all moment restrictions.
#' @param dir Search direction. \code{dir = 1} corresponds to looking for an
#' upper bound. \code{dir = -1} corresponds to looking for a lower bound.
#' @param EI.Mstep Function to compute expected improvements.
#' @param hyperparams List of hyperparameters.
#'
#' @returns Vector of starting values
get.starting.values <- function(theta.hash, dir, EI.Mstep, hyperparams) {

  # Extract necessary hyperparameters for readibility
  nbr.start.vals <- hyperparams[["nbr.start.vals"]]

  # Draw initial set of starting values
  draws.init <- draw.sv.init(theta.hash, dir, hyperparams, EI.Mstep)

  # Append theta.hash and modify them with the function MSpoints.R
  draws.init <- rbind(draws.init, c(theta.hash, 0))
  start.vals <- MSpoint(draws.init)

  # Only keep unique starting values
  start.vals <- start.vals[which.unique(matrix(start.vals, ncol = 1))]

  # Compute the expected improvement for each starting value
  EI.start.vals <- rep(0, length(start.vals))
  for (i in 1:length(start.vals)) {
    EI.start.vals[i] <- EI.Mstep(start.vals[i])
  }

  # Keep top nbr.start.vals starting values with highest expected improvement
  idxs.largest <- tail(order(EI.start.vals), nbr.start.vals)
  start.vals <- start.vals[idxs.largest]

  # Also include theta.hash and a slight perturbation of it
  theta.eps <- theta.hash + 1e-4
  start.vals <- c(start.vals, theta.hash, theta.eps)

  # Return the results
  start.vals
}

#' @title Optimize the expected improvement
#'
#' @description This function finds the point for which the expected improvement
#' is optimal, based on a given set of starting values. (M_step.R)
#'
#' @param start.vals Starting values for optimization.
#' @param EI.Mstep Function to compute expected improvements.
#' @param hyperparams List of hyperparameters.
#'
#' @returns Maximum of the expected imrpovement function.
#'
#' @import stats
#' @importFrom utils tail
#'
do.optimization.Mstep <- function(start.vals, EI.Mstep, hyperparams) {

  # Bounds of the parameter space for theta
  theta.lb <- hyperparams[["theta.lb"]]
  theta.ub <- hyperparams[["theta.ub"]]

  # Number of optimal values to return
  nbr.opt.EI <- hyperparams[["nbr.opt.EI"]]

  # Initialize object that will store the results of the optimization
  opt.EI <- rep(0, length(start.vals))
  opt.theta <- rep(0, length(start.vals))

  # Starting from each starting value, optimize the expected improvement
  for (i in 1:length(start.vals)) {

    # Set the starting value of this iteration
    start.val <- start.vals[i]

    # Maximize the expected improvement function
    opt.out <- optim(start.val, EI.Mstep, method = "L-BFGS-B", lower = theta.lb,
                     upper = theta.ub, control = list(fnscale = -1))
    opt.EI[i] <- opt.out$value
    opt.theta[i] <- opt.out$par
  }

  # Keep the top nbr.top.EI unique points with the highest expected improvement
  idxs.unique <- which.unique(matrix(opt.theta, ncol = 1))
  opt.EI <- opt.EI[idxs.unique]
  opt.theta <- opt.theta[idxs.unique]

  idxs.keep <- tail(order(opt.EI), nbr.opt.EI)
  opt.EI <- opt.EI[idxs.keep]
  opt.theta <- opt.theta[idxs.keep]

  # Return the results
  cbind(opt.theta, opt.EI)
}

#' @title Get extra evaluation points for E-step
#'
#' @description Function used to obtain extra theta values to be supplied to the
#' E-step in the next iteration (M_step.R). Note: this function should be
#' changed when implementing the sample space contractions (see comment made in
#' documentation of \code{M_step}).
#'
#' @param dir Search direction. \code{dir = 1} corresponds to looking for an
#' upper bound. \code{dir = -1} corresponds to looking for a lower bound.
#' @param theta.hash Tentative optimal value for theta, i.e., the largest or
#' smallest feasible value for theta (if dir = 1 or dir = -1, respectively). A
#' 'feasible value' is one that satisfies all moment restrictions.
#' @param maxviol.hash Violation curve evaluated at  \code{theta.hash}.
#' @param hyperparams List of hyperparameters.
#'
#' @returns Points to evaluate in E-step.
#'
#' @import stats
#'
get.extra.Estep.points <- function(dir, theta.hash, maxviol.hash,
                                   hyperparams) {

  # Extract the necessary information
  theta.lb <- hyperparams[["theta.lb"]]
  theta.ub <- hyperparams[["theta.ub"]]
  nbr.extra <- hyperparams[["nbr.extra"]]
  EAM_thetadistort <- hyperparams[["min.improvement"]]

  # Randomly sampled points
  rsp <- runif(nbr.extra, min = theta.lb, max = theta.ub)

  # Initialize object that will store perturbed theta.hash values
  theta_eps <- rep(0, 3)

  # Slight perturbations of theta.hash
  delta1 <- abs(maxviol.hash)
  delta2 <- EAM_thetadistort
  delta3 <- 10*EAM_thetadistort
  if (dir == 1) {
    theta_eps[1] <- min(theta.hash + delta1, theta.ub)
    theta_eps[2] <- min(theta.hash + delta2, theta.ub)
    theta_eps[3] <- min(theta.hash + delta3, theta.ub)
  } else if (dir == -1) {
    theta_eps[1] <- max(theta.hash - delta1, theta.lb)
    theta_eps[2] <- max(theta.hash - delta2, theta.lb)
    theta_eps[3] <- max(theta.hash - delta3, theta.lb)
  }

  # Return results
  c(rsp, theta_eps)
}

#### Search strategies: steps in EAM algorithm ####

#' @title Method for finding initial points of the EAM algorithm
#'
#' @description Also called the 'initialization' step in KMS19, this method
#' tries to find at least one initial feasible point, which is required to run
#' the EAM algorithm.
#' ToDo: Investigate whether the feasible point search of Bei (2024) is better.
#' If so, implement it.
#'
#' @param test.fun Function that takes a parameter vector as a first argument
#' and returns the test statistic, as well as the critical value.
#' @param hyperparams List of hyperparameters.
#' @param verbose Verbosity parameter.
#' @param picturose Picturosity flag. If \code{TRUE}, a plot illustrating the
#' workings of the algorithm will updated during runtime. Default is
#' \code{picturose = FALSE}.
#' @param parallel Flag for whether or not parallel computing should be used.
#' Default is \code{parallel = FALSE}.
#'
#' @returns Results of the initial feasible point search.
#'
#' @importFrom foreach foreach
#'
#' @references Kaido et al. (2019). Confidence intervals for projections of
#' partially identified parameters. Econometrica. 87(4):1397-1432.
#'
feasible_point_search <- function(test.fun, hyperparams, verbose,
                                  picturose = FALSE, parallel = FALSE) {

  #### Precondition checks ####

  # For ease of readability
  lb.theta <- hyperparams[["theta.lb"]]
  ub.theta <- hyperparams[["theta.ub"]]
  min.eval <- hyperparams[["min.eval"]]
  max.eval <- hyperparams[["max.eval"]]

  # Bounds on theta should be finite
  if ((abs(lb.theta) == Inf) | (abs(ub.theta) == Inf)) {
    stop("Bound on theta must be finite.")
  }

  # Minimum number of points to evaluate should be strictly positive
  if (min.eval <= 0) {
    stop("Minimum number of points to evaluate must be strictly positive.")
  }

  # Minimum number of evaluations should be smaller than maximum number of
  # evaluations
  if (min.eval > max.eval) {
    stop("min.eval should be smaller than max.eval.")
  }

  # If parallel computing is selected, a cluster should be initialized
  if (parallel) {
    if (!("clust" %in% names(hyperparams))) {
      stop(paste0("When parallel = TRUE, a parallel cluster should be ",
      "initialized and added to list of hyperparameters"))
    } else {
      clust <- hyperparams$clust
    }
  }

  #### Initial points to evaluate. If feasible point is found, return ####

  # Define initial set of points to evaluate
  pte <- seq(lb.theta, ub.theta, length.out = min.eval)

  # If plot of estimation procedure should be drawn, do so
  if (picturose) {
    plot_addpte(pte, col = "white")
  }

  # Initialize object that will store the test statistics and critical value for
  # each point in this initial set of points to evaluate.
  evaluations <- matrix(nrow = min.eval, ncol = 3)
  colnames(evaluations) <- c("theta", "t.stat", "crit.val")

  # For each point, obtain the test statistic and critical value.
  if (parallel) { # Using parallel computing

    # Update the user. Add to plot.
    if (verbose >= 3) {
      message("Evaluating initial points in parallel")
    }

    # Compute batch size and batch indices
    batch.size <- length(clust)
    n.batchs <- ceiling(length(pte)/batch.size)
    for (batch.idx in 1:n.batchs) {

      # Compute indices of pte vector to handle in this batch
      i.start <- (batch.idx - 1)*batch.size + 1
      i.end <- min(batch.idx * batch.size, length(pte))

      # If required, add points to plot
      if (picturose) {
        plot_addpte(pte[i.start:i.end])
      }

      # Evaluate batch of points.
      suppressWarnings({evaluations[i.start:i.end, ] <-
        foreach(i = i.start:i.end, .combine = 'rbind',
                .export = c("par.space", "hp", "c", "inst.func.evals", "t",
                            "data", "test.fun", "options")
        ) %dopar% {

          # Load all necessary packages
          source("simulationFunctions.R")

          # Select theta of this iteration
          theta <- pte[i]

          # Run the test
          test.out <- test.fun(theta)

          # Return the results
          c(theta, test.out[["t.stat"]], test.out[["crit.val"]])
        }})

      # If required, add points to plot
      if (picturose) {
        plot_addpte.eval(evaluations[i.start:i.end, , drop = FALSE])
      }
    }

  } else { # Using sequential computing

    for (i in 1:min.eval) {
      if (verbose >= 3) {
        message(sprintf("Checking initial points (%s / %s)", i, min.eval))
      }

      # Select theta of this iteration
      theta <- pte[i]

      # If necessary, plot it
      if (picturose) {
        plot_addpte(theta)
      }

      # Run the test
      test.out <- test.fun(theta)

      # Store the results
      evaluations[i,] <- c(theta, test.out[["t.stat"]], test.out[["crit.val"]])

      # If plot of estimation procedure should be drawn, do so.
      if (picturose) {
        plot_addpte.eval(evaluations[i, , drop = FALSE])
      }
    }
  }

  # If there is at least one feasible point in the set of evaluated points,
  # return the set of evaluated points.
  if (any(evaluations[, "t.stat"] <= evaluations[, "crit.val"])) {
    return(list(evaluations = evaluations))
  }

  #### If no feasible point was found, continue search ####

  if (parallel) { # Using parallel computing

    # Set some parameters
    batch.size <- length(clust)
    eval.nbr <- min.eval
    stop <- FALSE

    # Evaluate (in batch) additional points until stopping criterion is reached
    while(eval.nbr < max.eval & !stop) {

      # Obtain next batch of points to evaluate
      evaluations.dummy <- evaluations
      pte <- c()
      pte.idxafter <- c()
      for (i in 1:batch.size) {
        gnp.out <- get.next.point(evaluations.dummy, lb.theta, ub.theta)
        theta.next <- gnp.out[["theta.next"]]
        idx.after <- gnp.out[["idx.after"]]
        row <- c(theta.next, 1, 0)
        evaluations.dummy <- insert.row(evaluations.dummy, row, idx.after)
        pte <- c(pte, theta.next)
        pte.idxafter <- c(pte.idxafter, idx.after)
      }

      # Evaluate batch
      evaluations.add <-
        foreach(i = 1:batch.size, .combine = 'rbind') %dopar% {

          # Load all necessary packages
          source("simulationFunctions.R")

          # Select theta of this iteration
          theta <- pte[i]

          # Run the test
          test.out <- test.fun(theta)

          # Return the results
          c(theta, test.out[["t.stat"]], test.out[["crit.val"]])
        }

      # Store the results
      for (i in 1:batch.size) {
        evaluations <- insert.row(evaluations, evaluations.add[i, ], pte.idxafter[i])
      }

      # Update stopping criteria
      stop <- any(evaluations[, "t.stat"] < evaluations[, "crit.val"])
      eval.nbr <- eval.nbr + batch.size
    }

  } else { # Using sequential computing

    # Initialization for while-loop
    gnp.out <- get.next.point(evaluations, lb.theta, ub.theta)
    theta.next <- gnp.out[["theta.next"]]
    idx.after <- gnp.out[["idx.after"]]
    stop <- gnp.out[["stop"]]
    eval.nbr <- min.eval + 1

    # If no feasible point has been found yet, test if the midpoints
    # between evaluated points are feasible. Continue this procedure until either
    # a feasible point is found, or the the stopping criterion is reached.
    while (!stop & (eval.nbr <= max.eval)) {

      # Update user
      if (verbose >= 3) {
        message(sprintf("Checking additional points (%s / %s)", eval.nbr, max.eval))
      }
      if (picturose) {
        plot_addpte(theta.next)
      }

      # Run the test for this point
      test.out <- test.fun(theta.next)

      # Add the results to the set of evaluated points
      row <- c(theta.next, test.out[["t.stat"]], test.out[["crit.val"]])
      evaluations <- insert.row(evaluations, row, idx.after)

      # Get next point
      gnp.out <- get.next.point(evaluations, lb.theta, ub.theta)
      theta.next <- gnp.out[["theta.next"]]
      idx.after <- gnp.out[["idx.after"]]
      stop <- gnp.out[["stop"]]

      # Increment number of evaluations
      eval.nbr <- eval.nbr + 1

      # Update plot
      if (picturose) {
        plot_addpte.eval(evaluations[nrow(evaluations), , drop = FALSE])
      }
    }

  }

  #### Return the results ####

  list(evaluations = evaluations)
}

#' @title E-step in the EAM algorithm as described in KMS19.
#'
#' @description This function performs the estimation step in the EAM algorithm.
#'
#' @param thetas Points at which to perform the E-step. Usually the result of
#' the M-step.
#' @param test.fun Function returning the test statistic, as well as the critical
#' value.
#' @param dir Direction in which to optimize. For finding upper bounds, set
#' \code{dir = 1}, for finding lower bounds, set \code{dir = -1}.
#' @param evaluations Matrix containing each point that was already evaluated,
#' alongside the corresponding test statistic and critical value, as its rows.
#' @param verbose Verbosity parameter.
#'
#' @returns Results of the E-step.
#'
E_step <- function(thetas, test.fun, dir, evaluations, verbose) {

  # Loop over all values of theta to be checked
  for (i in 1:length(thetas)) {

    # Get the theta value of this iteration
    theta <- thetas[i]

    # Check whether the evaluation of the current theta value was already
    # carried out and can hence be skipped
    if (any(abs(theta - evaluations[, "theta"]) < 1e-10)) {
      if (verbose >= 3) {
        message(sprintf(
          "\t Evaluating point %s out of %s... (skipped)", i, length(thetas)
        ))
      }
      next
    }

    if (verbose >= 3) {
      message(sprintf("\t Evaluating point %s out of %s...", i, length(thetas)))
    }

    # Obtain the test statistic and critical value of the given point
    test.out <- test.fun(theta)
    t.stat <- test.out[["t.stat"]]
    crit.val <- test.out[["crit.val"]]

    # Append the point to the set of evaluated points
    evaluations <- rbind(evaluations, c(theta, t.stat, crit.val))
    evaluations <- evaluations[order(evaluations[, 1]),]
  }

  # Indices of feasible points
  feas.idxs <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])

  # Tentative optimal value
  if (dir == 1) {
    theta.hash <- max(evaluations[feas.idxs, "theta"])
  } else if (dir == -1) {
    theta.hash <- min(evaluations[feas.idxs, "theta"])
  }

  # Return the results
  list(evaluations = evaluations, theta.hash = theta.hash)
}

#' @title A-step in the EAM algorithm described in KMS19
#'
#' @description This function performs the approximation step in the EAM
#' algorithm. More specifically, it fits a Gaussian-process regression model
#' (Kriging) to the evaluated data points \eqn{(\theta, c(\theta))}.
#'
#' @param evaluations Matrix containing each point that was already evaluated,
#' alongside the corresponding test statistic and critical value, as its rows.
#' @param verbose Verosity parameter.
#'
#' @returns Results of the A-step.
#'
#' @importFrom graphics abline lines
#' @importFrom utils install.packages
#' @seealso Package \pkg{SPOT} on \url{https://CRAN.R-project.org/package=SPOT}.
#'
A_step <- function(evaluations, verbose = 0) {

  # Make sure package SPOT is installed and loaded
  while (!requireNamespace("SPOT", quietly = TRUE)) {
    ans <- readline("Package SPOT is required for EAM algorithm but not available. Would you like to install it? (Y/N)")
    if (regexpr(ans, 'y', ignore.case = TRUE) == 1) {
      install.packages("SPOT")
    } else if (regexpr(ans, 'n', ignore.case = TRUE) == 1) {
      stop("Unable to run EAM algorithm. Please specify a different root finding algorithm.")
    } else {
      message("Invalid answer. Please only answer with 'y' (yes) or 'n' (no). Returning to original input line...")
    }
  }

  # Matrix of theta values
  theta.mat <- matrix(evaluations[, "theta"], nrow = nrow(evaluations))

  # Make sure design points are not too close together, as it will cause
  # unwanted behaviour of the kriging model.
  idxs.unique <- which.unique(theta.mat, tol = 1/nrow(theta.mat))
  theta.mat <- theta.mat[idxs.unique, , drop = FALSE]

  # Vector of violations (test statistic - critical value)
  violations.vct <- evaluations[idxs.unique, "t.stat"] - evaluations[idxs.unique, "crit.val"]

  # Control parameters
  control = list(regr = SPOT::regpoly0, corr = SPOT::corrkriging, target = c("y", "s"))

  # Fit the Kriging model
  fit.krige <- SPOT::buildKrigingDACE(theta.mat, violations.vct, control = control)

  # If asked, plot the Kriging model
  if (verbose >= 3) {
    x.vals <- seq(min(evaluations[, "theta"]), max(evaluations[, "theta"]),
                  length.out = 500)
    predictions <- predict(fit.krige, matrix(x.vals, nrow = length(x.vals)))
    y.vals <- predictions$y
    sd.vals <- predictions$s

    plot(x.vals, y.vals, type = 'l', xlab = "theta", ylab = "predicted violation",
         main = "Kriging model")
    lines(x.vals, y.vals + 2 * sd.vals, type = 'l', lty = 2, col = "red")
    lines(x.vals, y.vals - 2 * sd.vals, type = 'l', lty = 2, col = "red")
    abline(h = 0, col = "blue", lty = 2)
  }

  # Return the Kriging model
  fit.krige
}

#' @title M-step in the EAM algorithm described in KMS19.
#'
#' @description  This function performs the maximization step in the EAM
#' algorithm. More specifically, it maximizes the expected improvement.
#' ToDo: implement sample space contractions (see comment made in documentation
#' of \code{draw.sv.init}).
#'
#' @param dir Direction to search in. \code{dir = 1} corresponds to finding the
#' upper bound of the confidence interval. \code{dir = -1} corresponds to
#' finding the lower bound.
#' @param evaluations Matrix containing each point that was already evaluated,
#' alongside the corresponding test statistic and critical value, as its rows.
#' @param theta.hash Tentative best value of theta. Obtained from the E-step.
#' @param fit.krige Kriging model obtained from the A-step.
#' @param test.fun The test function to be inverted in order to obtain the
#' identified set.
#' @param c Projection vector.
#' @param par.space Bounds of the parameter space.
#' @param hyperparams Parameters used in obtaining initial values
#' for the maximization algorithm. If \code{NULL}, default values are used.
#' Default is \code{hyperparams = NULL}.
#' @param verbose Verbosity parameter.
#'
#' @importFrom graphics abline
#'
M_step <- function(dir, evaluations, theta.hash, fit.krige, test.fun, c,
                   par.space, hyperparams, verbose) {

  #### Set some useful variables ####

  # Compute tentative optimal value for theta if not supplied
  if (is.null(theta.hash)) {

    # Indices of feasible points
    feas.idxs <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])

    # Find current best value
    if (dir == 1) {
      theta.hash <- max(evaluations[feas.idxs, "theta"])
    } else if (dir == -1) {
      theta.hash <- min(evaluations[feas.idxs, "theta"])
    }
  }

  # Current largest value of feasible theta
  theta.hash.idx <- which(evaluations[, "theta"] == theta.hash)[1]
  maxviol.hash <- evaluations[theta.hash.idx, "t.stat"] - evaluations[theta.hash.idx, "crit.val"]

  # Define expected improvement function wrt theta.hash
  EI.Mstep <- function(theta) {EI(theta, test.fun, fit.krige, theta.hash, dir)}

  #### Main code for M-step ####

  # Obtain starting values for the maximization of the EI
  start.vals <- get.starting.values(theta.hash, dir, EI.Mstep, hyperparams)

  # Compute the optimal expected improvement and corresponding theta-value.
  opt.res <- do.optimization.Mstep(start.vals, EI.Mstep, hyperparams)

  # Add a randomly drawn theta s.t. dir * theta > dir * theta.hash to the set
  # of points to evaluate in the next E-step.
  extra.points <- get.extra.Estep.points(dir, theta.hash, maxviol.hash,
                                         hyperparams)

  # If asked, plot the expected improvement function, theta.hash and the most
  # promising theta value.
  if (verbose >= 3) {
    x.vals <- seq(min(evaluations[, "theta"]), max(evaluations[, "theta"]),
                  length.out = 500)
    y.vals <- Vectorize(EI.Mstep)(x.vals)

    plot(x.vals, y.vals, type = 'l', xlab = "theta", ylab = "EI",
         main = "Expected improvement function M-step")
    abline(v = theta.hash, col = "red")
    abline(v = opt.res[1], col = "green")
  }

  # Return the results
  list(opt.res = opt.res, extra.points = extra.points)
}

#' @title Check convergence of the EAM algorithm.
#'
#' @description This function checks the convergence of the EAM algorithm.
#' ToDo: Get rid of hard coding stop of algorithm when no more improvement of
#' theta values (maybe related to parameter space contraction, since the problem
#' is that the given points to check in the E-step of the following iteration
#' can always be the same and always be rejected (except of course for the
#' randomly chosen one), while the most promising theta value continues to be
#' the same, infeasible value. In this way, it is possible that
#' theta.hash - mp.theta.next at some point will never decrease).
#'
#' @param opt.val.prev Previous optimal theta value.
#' @param evaluations Matrix of violation curve evaluations.
#' @param mp.theta.next Most promising value of theta for which to run the
#' E-step in the following iteration
#' @param iter.nbr Number of iterations of the EAM algorithm run so far.
#' @param dir Search direction.
#' @param hyperparams List of hyperparameters used in the EAM algorithm.
#' @param verbose Verbosity parameter.
#'
#' @returns Boolean value whether or not algorithm has converged.
#'
EAM.converged <- function(opt.val.prev, evaluations, mp.theta.next, iter.nbr,
                          dir, hyperparams, verbose) {

  #### Extract necessary information from set of hyperparameters ####

  # Minimum amount that this step should have improved upon the previous step
  min.improvement <- hyperparams[["min.improvement"]]

  # Minimum amount of improvement that can be the result of running the next
  # step.
  min.possible.improvement <- hyperparams[["min.possible.improvement"]]

  # Minimum number of EAM iterations that should be run.
  EAM.min.iter <- hyperparams[["EAM.min.iter"]]

  # Maximum number of EAM iterations that should be run.
  max.iter <- hyperparams[["max.iter"]]

  #### Check all convergence criteria ####

  # Indices of feasible points
  feas.idxs <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])

  # Best and second best feasible point
  if (dir == 1) {
    opt.val <- sort(evaluations[feas.idxs, "theta"], decreasing = TRUE)[1]
  } else if (dir == -1) {
    opt.val <- sort(evaluations[feas.idxs, "theta"], decreasing = FALSE)[1]
  }

  # If this is run in the first iteration, opt.val.prev = NULL. In this case,
  # set opt.val.prev = -Inf
  if (is.null(opt.val.prev)) {
    opt.val.prev <- -Inf
  }

  conv1 <- (dir * (opt.val - opt.val.prev) < min.improvement)
  conv2 <- (dir * (mp.theta.next - opt.val) < min.possible.improvement)
  conv3 <- (iter.nbr > EAM.min.iter)

  # Determine whether or not to stop the algorithm (stopping criterion based on
  # implementation of Bei, 2024).
  stop <- (conv1 & conv2 & conv3) | (iter.nbr > max.iter)

  # If there is no more improvement in theta value, stop the algorithm. This is
  # valid, since each iteration tries to improve the current best theta value by
  # 'min.improvement'.
  #
  # !!! See also [ToDo] !!!
  if (abs(opt.val - opt.val.prev) < 1e-12) {
    stop <- TRUE
  }

  # If necessary, print convergence information to the console
  if (verbose >= 3) {
    message(
      "_______________________________________________________________________________")
    message(sprintf("Iteration %s:\n", iter.nbr))
    message(
      "Impr. wrt previous | Possible improvement | min.iter reached | max.iter reached"
    )
    message(
      "-------------------|----------------------|------------------|-----------------"
    )
    message(sprintf(
      "%18f | %20f | %16s | %15s", dir * (opt.val - opt.val.prev),
      dir * (mp.theta.next - opt.val), ifelse(conv3, "TRUE", "FALSE"),
      ifelse(iter.nbr > max.iter, "TRUE", "FALSE")
    ))
    message(
      "_______________________________________________________________________________")
  } else if (verbose %in% c(1, 2)) {
    if (iter.nbr == 1) {
      message(
        "______________________________________________________________________________________"
      )
      message(
        "iter | Impr. wrt previous | Possible improvement | min.iter reached | max.iter reached"
      )
      message(
        "-----|--------------------|----------------------|------------------|-----------------"
      )
    }
    message(sprintf(
      "%4d | %18f | %20f | %16s | %15s", iter.nbr, dir * (opt.val - opt.val.prev),
      dir * (mp.theta.next - opt.val), ifelse(conv3, "TRUE", "FALSE"),
      ifelse(iter.nbr > max.iter, "TRUE", "FALSE")
    ))
  }

  # Return the result
  stop
}

#### Search strategies: steps in grid search ####

#' @title Rudimentary, bidirectional 1D grid search algorithm.
#'
#' @description
#' This function implements a rudimentary, bidirectional search algorithm. It
#' works by expanding a grid with given step.size in both directions, starting
#' from an initial feasible point.
#'
#' @param test.results Matrix containing the evaluations of the test statistic
#' and critical value.
#' @param max.iter Maximum number of iterations.
#' @param step.size Step size based on which the grid is constructed.
#'
#' @returns The next point to evaluate in the grid search.
#'
gs.algo.bidir <- function(test.results, max.iter, step.size) {

  # Define some useful variables
  n.test <- nrow(test.results)

  # Point around which the grid is centered
  center.point <- test.results[1, 1]

  # Initialize some variables
  stop <- FALSE

  #### Some edge cases ####

  # If test.results only contains the initial point, confirm that the test
  # did not lead to a rejection and if so, return a new r value to the right of
  # the original one. Else stop the algorithm
  if (n.test == 1) {
    if (test.results[n.test, 2] <= test.results[n.test, 3]) {
      r <- center.point + step.size
      return(list(r = r, stop = stop))
    } else {
      return(list(r = NULL, stop = TRUE))
    }
  }

  # If test.results contains the initial point and one initial point, return an
  # r evaluation point in the other direction of the initial point.
  if (n.test == 2) {
    r <- center.point - sign(test.results[n.test, 1] - center.point) * step.size
    return(list(r = r, stop = stop))
  }

  # If maximum number of grid evaluations is reached, stop the algorithm
  if (n.test >= max.iter) {
    return(list(r = NULL, stop = TRUE))
  }

  #### Main logic ####

  # If the test evaluations for the last two r values where at opposite sides of
  # the center point...
  if ((test.results[n.test, 1] - center.point) * (test.results[n.test - 1, 1] - center.point) < 0) {

    # If the second to last r value led to a non-rejection...
    if (test.results[n.test - 1, 2] <= test.results[n.test - 1, 3]) {

      # Return an r value one step further than this second to last r value
      r <- test.results[n.test - 1, 1] + sign(test.results[n.test - 1, 1]) * step.size
    }

    # If the second to last r value led to a rejection...
    else {

      # Check whether the last r value led to a rejection. If it didn't lead to
      # a rejection...
      if (test.results[n.test, 2] <= test.results[n.test, 3]) {

        # Return an r value one step further than the last r value
        r <- test.results[n.test, 1] + sign(test.results[n.test, 1]) * step.size
      }

      # If it did lead to a rejection...
      else {

        # Stop the algorithm
        r <- NULL
        stop <- TRUE

      }
    }
  }

  # If the test evaluation for the last two r values was at the same side of the
  # center point...
  else {

    # Test whether the last r value evaluation led to a rejection. If it didn't...
    if (test.results[n.test, 2] <= test.results[n.test, 3]) {

      # Return an r value one step further than the last r value
      r <- test.results[n.test, 1] + sign(test.results[n.test, 1]) * step.size
    }

    # If if did...
    else {

      # Stop the algorithm
      r <- NULL
      stop <- TRUE
    }
  }

  # Return results
  return(list(r = r, stop = stop))
}

#' @title Return the next point to evaluate when doing regular grid search
#'
#' @description This function implements a unidirectional grid search, that
#' works by expanding a grid starting from a given feasible point in the
#' given direction.
#'
#' @param evaluations Matrix of evaluated test statistics and critical values.
#' @param dir Search direction.
#' @param iter.nbr Iteration number.
#' @param hp List of hyperparameters.
#'
#' @returns Next point to evaluate in the search algorithm.
gs.regular <- function(evaluations, dir, iter.nbr, hp) {

  # Extract the necessary hyperparameters
  step.size <- hp[["step.size"]]
  max.iter <- hp[["max.iter"]]

  # Find current largest (smallest) feasible value, called theta.hash
  feas.idx <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])
  theta.hash <- dir * max(dir * evaluations[feas.idx, "theta"])

  # If theta.hash is larger (smaller) than the upper (lower) bound, stop the
  # algorithm.
  if (((dir == 1) & (theta.hash > hp$theta.ub)) |
      ((dir == -1) & (theta.hash < hp$theta.lb))) {
    theta.to.eval <- theta.hash + dir
    return(list(theta.to.eval = theta.to.eval, stop = TRUE))
  }

  # Get next theta value
  theta.to.eval <- theta.hash + dir * step.size

  # Check whether theta.to.eval has already been evaluated and if so, whether or
  # not it was feasible.
  stop <- FALSE
  if (theta.to.eval %in% evaluations[, "theta"]) {
    idx.tte <- which(evaluations[, "theta"] == theta.to.eval)
    stop <- evaluations[idx.tte, "t.stat"] > evaluations[idx.tte, "crit.val"]
  }
  if (iter.nbr > max.iter) {
    stop <- TRUE
  }

  # Return the results
  list(theta.to.eval = theta.to.eval, stop = stop)
}

#' @title Return the next point to evaluate when doing binary search
#'
#' @description This function implements the binary search algorithm, that
#' starts from a given feasible point and looks in the given direction for the
#' root of the violation curve.
#'
#' @param evaluations Matrix of evaluated test statistics and critical values.
#' @param dir Search direction.
#' @param iter.nbr Iteration number.
#' @param hp List of hyperparameters.
#'
#' @returns The next point to evaluate.
#'
gs.binary <- function(evaluations, dir, iter.nbr, hp) {

  # Extract the necessary hyperparameters
  bin.search.tol <- hp[["bin.search.tol"]]
  max.iter <- hp[["max.iter"]]

  # Check if maximum number of iterations has been reached. If so, stop the
  # algorithm
  if (iter.nbr > max.iter) {
    return(list(theta.to.eval = NULL, stop = TRUE))
  }

  # Find current largest (smallest) feasible value, called theta.hash
  feas.idx <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])
  theta.hash <- dir * max(dir * evaluations[feas.idx, "theta"])

  # If theta.hash is larger (smaller) than or equal to the upper (lower) bound,
  # stop the algorithm.
  if (((dir == 1) & (theta.hash >= hp$theta.ub)) |
      ((dir == -1) & (theta.hash <= hp$theta.lb))) {
    return(list(theta.to.eval = NULL, stop = TRUE))
  }

  # Matrix of all infeasible points
  infeas.idxs <- which(evaluations[, "t.stat"] > evaluations[, "crit.val"])
  evals.infeas <- matrix(nrow = 0, ncol = 3)
  colnames(evals.infeas) <- colnames(evaluations)
  for (idx in infeas.idxs) {
    evals.infeas <- rbind(evals.infeas, evaluations[idx, ])
  }

  # Indices of points in infeasible point set that are larger (smaller) than
  # theta.hash
  idxs <- which(dir * evals.infeas[, "theta"] > dir * theta.hash)

  # If no infeasible point which is larger (smaller) than theta.hash is found
  # yet, search for one.
  #
  # NOTE: when both theta.lb and theta.ub are checked in an initial stage, this
  #       code block will never be ran. Indeed, either theta.lb (theta.ub) is
  #       feasible, in which case the root finding algorithm would already
  #       have been stopped (checked above). If it is infeasible, then an
  #       infeasible point that is smaller (larger) than theta.hash is known.
  iter.nbr <- 1
  if (length(idxs) == 0) {

    # Obtain largest (smallest) and second largest (smallest) theta values
    out.sort <- dir * sort(dir * evaluations[, "theta"], decreasing = TRUE)[1:2]
    theta.largest <- out.sort[1]
    theta.second <- out.sort[2]

    # Theta value to evaluate
    dist <- 2 * max(abs(theta.largest - theta.second), 1)
    theta.to.eval <- dir * (dir * theta.largest + dist)

    # Return theta
    return(list(theta.to.eval = theta.to.eval, stop = FALSE))
  }

  # Else, determine theta.tilde
  theta.tilde <- dir * min(dir * evals.infeas[idxs, "theta"])

  # Midpoint between theta.tilde and theta.hash
  theta.to.eval <- (theta.hash + theta.tilde)/2

  # Stopping criterion
  stop <- FALSE
  if (abs(theta.to.eval - theta.hash) < bin.search.tol) {
    stop <- TRUE
  }

  # Return the results
  list(theta.to.eval = theta.to.eval, stop = stop)
}

#' @title Return the next point to evaluate when doing interpolation search
#'
#' @description This function implements the interpolation search algorithm,
#' that starts from a given feasible point and looks in the given direction for
#' the root of the violation curve.
#'
#' @param evaluations Matrix of evaluated test statistics and critical values.
#' @param dir Search direction.
#' @param iter.nbr Iteration number.
#' @param hp List of hyperparameters.
#'
#' @returns The next point to evaluate.
#'
gs.interpolation <- function(evaluations, dir, iter.nbr, hp) {

  # Extract the necessary hyperparameters
  int.search.tol <- hp[["bin.search.tol"]]
  max.iter <- hp[["max.iter"]]

  # Check if maximum number of iterations has been reached. If so, stop the
  # algorithm
  if (iter.nbr > max.iter) {
    return(list(theta.to.eval = NULL, stop = TRUE))
  }

  # Find current largest (smallest) feasible value, called theta.hash
  feas.idx <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])
  theta.hash <- dir * max(dir * evaluations[feas.idx, "theta"])

  # If theta.hash is larger (smaller) than or equal to the upper (lower) bound,
  # stop the algorithm.
  if (((dir == 1) & (theta.hash >= hp$theta.ub)) |
      ((dir == -1) & (theta.hash <= hp$theta.lb))) {
    return(list(theta.to.eval = NULL, stop = TRUE))
  }

  # Matrix of all infeasible points
  infeas.idxs <- which(evaluations[, "t.stat"] > evaluations[, "crit.val"])
  evals.infeas <- matrix(nrow = 0, ncol = 3)
  colnames(evals.infeas) <- colnames(evaluations)
  for (idx in infeas.idxs) {
    evals.infeas <- rbind(evals.infeas, evaluations[idx, ])
  }

  # Indices of points in infeasible point set that are larger (smaller) than
  # theta.hash
  idxs <- which(dir * evals.infeas[, "theta"] > dir * theta.hash)

  # If no infeasible point which is larger (smaller) than theta.hash is found
  # yet, search for one.
  #
  # NOTE: when both theta.lb and theta.ub are checked in an initial stage, this
  #       code block will never be ran. Indeed, either theta.lb (theta.ub) is
  #       feasible, in which case the root finding algorithm would already
  #       have been stopped (checked above). If it is infeasible, then an
  #       infeasible point that is smaller (larger) than theta.hash is known.
  iter.nbr <- 1
  if (length(idxs) == 0) {

    # Obtain largest (smallest) and second largest (smallest) theta values
    out.sort <- dir * sort(dir * evaluations[, "theta"], decreasing = TRUE)[1:2]
    theta.largest <- out.sort[1]
    theta.second <- out.sort[2]

    # Theta value to evaluate
    dist <- 2 * max(abs(theta.largest - theta.second), 1)
    theta.to.eval <- dir * (dir * theta.largest + dist)

    # Return theta
    return(list(theta.to.eval = theta.to.eval, stop = FALSE))
  }

  # Else, determine theta.tilde
  theta.tilde <- dir * min(dir * evals.infeas[idxs, "theta"])
  idx.theta.tilde <- which(evals.infeas[, "theta"] == theta.tilde)

  # Obtain the values of the violation curve at both points
  theta.tilde.viol <- evals.infeas[idx.theta.tilde, "t.stat"] - evals.infeas[idx.theta.tilde, "crit.val"]
  idx.theta.hash <- which(evaluations[, "theta"] == theta.hash)
  theta.hash.viol <- evaluations[idx.theta.hash, "t.stat"] - evaluations[idx.theta.hash, "crit.val"]
  theta.to.eval <- theta.hash - theta.hash.viol * (theta.tilde - theta.hash)/(theta.tilde.viol - theta.hash.viol)

  # Stopping criterion
  stop <- FALSE
  if (abs(theta.to.eval - theta.hash) < int.search.tol) {
    stop <- TRUE
  }

  # Return the results
  list(theta.to.eval = theta.to.eval, stop = stop)
}

#### Search strategies: main functions EAM ####

#' @title Set default hyperparameters for EAM algorithm
#'
#' @description This function returns a list with the (default) hyperparameters
#' used in the EAM algorithm
#'
#' @param options A list of user-specified values for (some of) the
#' hyperparameters. These hyperparameters can include:
#' \describe{
#'  \item{min.dist/max.dist:}{The minimum/maximum distance of sampled points
#'  from the current best value for the coefficient of interest.}
#'  \item{min.eval/max.eval:}{The minimum/maximum number of points evaluated
#'  in the initial feasible point search.}
#'  \item{nbr.init.sample.points:}{The total number of drawn points required in
#'  the initial drawing process.}
#'  \item{nbr.init.unif:}{The total number of uniformly drawn points in the
#'  initial set of starting values.}
#'  \item{nbr.points.per.iter.init:}{Number of points sampled per iteration in
#'  the initial drawing process.}
#'  \item{nbr.start.vals:}{Number of starting values for which to run the
#'  optimization algorithm for the expected improvement.}
#'  \item{nbr.opt.EI:}{Number of optimal theta values found by the optimization
#'  algorithm to return.}
#'  \item{nbr.extra:}{Number of extra randomly drawn points to add to the set
#'  of optimal theta values (to be supplied to the next E-step).}
#'  \item{min.improvement:}{Minimum amount that the current best root of the
#'  violation curve should improve by wrt. the its previous value.}
#'  \item{min.possible.improvement:}{Minimum amount that the next iteration
#'  should be able to improve upon the current best value of the root.}
#'  \item{EAM.min.iter:}{Minimum amount of EAM iterations to run.}
#'  \item{max.iter:}{Maximum amount of EAM iterations to run.}
#' }
#'
#' @returns List of hyperparameters for the EAM algotithm.

set.EAM.hyperparameters <- function(options) {

  # Define the list of hyperparameters
  hyperparams <- list(

    # Minimum and maximum distance of sampled points from current best theta
    min.dist = 1e-4,
    max.dist = 1,

    # Minimum and maximum number of points to evaluate in initial feasible
    # search.
    min.eval = 10,
    max.eval = 100,

    # Total number of drawn points required in initial drawing process
    nbr.init.sample.points = 10,

    # Total number of uniformly drawn points in the initial set of starting
    # values
    nbr.init.unif = 20,

    # Number of points sampled per iteration in the initial drawing process
    nbr.points.per.iter.init = 4,

    # Number of starting values with which to run the optimization algorithm
    # for the expected improvement.
    nbr.start.vals = 80,

    # Number of optimal theta values found by the optimization algorithm to
    # return
    nbr.opt.EI = 1,

    # Number of extra randomly drawn points to add to the set of optimal
    # theta values (to be supplied to the next E-step).
    nbr.extra = 1,

    # [NOT USED] Distortion to be applied to theta.hash in order to obtain extra
    # evaluation points for the next E-step. In the code, 'min.improvement' is
    # used as the distortion value for theta.
    EAM_thetadistort = 0.005,

    # Minimum amount that this step should have improved upon the previous step
    min.improvement = 0.0001,

    # Minimum amount of improvement that can be the result of running the next
    # step.
    min.possible.improvement = 0.005,

    # Minimum number of EAM iterations that should be run
    EAM.min.iter = 4,

    # Maximum number of EAM iterations that should be run
    max.iter = 100
  )

  # Overwrite default values with user specified values if needed
  hyperparams[names(options)] <- options

  # Return the result
  hyperparams
}

#' @title Main function to run the EAM algorithm
#'
#' @description This function implements the EAM search strategy as described in
#' Kaido, Molinari and Stoye (2019). This is a generic function, requiring the
#' specification of a test function (\code{test.fun}), as well as the
#' specification of the parameter space (\code{hyperparams}).
#'
#' @param dir The direction of the test. \code{dir = 1} corresponds to finding
#' the upper bound of the identified set, \code{dir = -1} corresponds to finding
#' the lower bound.
#' @param test.fun The test function to be inverted in order to obtain the
#' identified set. It should take a scalar parameter as argument (i.e. the
#' specified value of a component of the full parameter vector) and return a
#' list with named elements \code{list(theta, t.stat, crit.val)}, where
#' \code{theta} is the scalar value that was tested, \code{t.stat} is the value
#' of the test statistic and \code{crit.val} is the critical value to be used in
#' determining whether to reject or not reject.
#' @param hyperparams A list of hyperparameters needed in order for this method
#' to run (see \code{set.EAM.hyperparameters.R}).
#' @param evaluations Matrix of already evaluated points, of which at least one
#' is feasible. When \code{evaluations = NULL} (default), the initial feasible
#' point search will be executed first.
#' @param time.run.duration Boolean value indicating whether to time each step
#' in the EAM algorithm. Requires \code{chronometer.R}. Default is
#' \code{time.run.duration = FALSE}.
#' @param verbose Boolean value indicating whether or not to print run time
#' updates to the console. Default is \code{verbose = FALSE}.
#' @param picturose Boolean value indicating whether or not to visualize the
#' identified set search. Default is \code{picturose = FALSE}.
#'
#' @returns List containing the evaluations of the test statistic and critical
#' values, convergence information, and run times.
#'
#' @references Kaido et al. (2019). Confidence intervals for projections of
#' partially identified parameters. Econometrica. 87(4):1397-1432.
#'
EAM <- function(dir, test.fun, hyperparams, evaluations = NULL,
                time.run.duration = FALSE, verbose = 0, picturose = FALSE) {

  # Extract parameter space
  par.space <- hyperparams$par.space

  # If the run times should be recorded, initialize a chronometer object
  if (time.run.duration) {
    source("chronometer.R")
    chronometer <- Chronometer$new()
    chronometer$start()
  } else {
    chronometer <- NULL
  }

  # Returns a set of points, of which at least one will be feasible
  if (!is.null(evaluations)) {
    fps.fail_flag <- FALSE
    if (all(evaluations[, "t.stat"] > evaluations[, "crit.val"])) {
      fps.fail_flag <- TRUE
    }
  } else {
    fps.out <- feasible_point_search(test.fun, hyperparams, verbose)
    evaluations <- fps.out[["evaluations"]]
    fps.fail_flag <- fps.out[["stop"]]
    if (time.run.duration) {
      chronometer$record.leg("finished.fps")
    }
  }

  # Get the tentative best value for theta
  feas.idxs <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])
  theta.hash <- dir * max(dir * evaluations[feas.idxs, "theta"])

  # Initialize variables used in the main loop
  iter.nbr <- 1               # The index of the current iteration
  converged <- fps.fail_flag  # Convergence flag
  ptc.Estep <- c()            # Vector of points to check in the E-step
  opt.theta.prev <- NULL      # Optimal theta value of previous iteration

  # Perform main loop
  if (verbose >= 3) {
    message("Entering main loop...")
  }
  while (!converged) {

    #### Evaluation step ####

    # Run the test on the theta values returned in the previous iteration.
    # (This step is skipped during the first iteration)
    if (iter.nbr > 1) {
      if (verbose >= 3) {
        message("Doing E-step...")
      }
      E_step.out <- E_step(ptc.Estep, test.fun, dir, evaluations, verbose)
      evaluations <- E_step.out[["evaluations"]]
      theta.hash <- E_step.out[["theta.hash"]]
    }

    if (time.run.duration) {
      chronometer$record.leg("finished.Estep")
    }

    #### Approximation step ####

    if (verbose >= 3) {
      message("Doing A-step...")
    }

    # Fit a Kriging model using the points (theta, t.stat(theta) - c(theta)).
    fit.krige <- A_step(evaluations, verbose)

    if (time.run.duration) {
      chronometer$record.leg("finished.Astep")
    }

    #### Maximization step ####

    if (verbose >= 3) {
      message("Doing M-step...")
    }

    # Find promising value(s) of theta to be checked in the next iteration of
    # this algorithm. To counter greediness of the search, also include some
    # randomly selected points.
    M_step.out <- M_step(dir, evaluations, theta.hash, fit.krige, test.fun, c,
                         par.space, hyperparams, verbose)
    opt.res <- M_step.out[["opt.res"]]
    extra.points <- M_step.out[["extra.points"]]

    if (time.run.duration) {
      chronometer$record.leg("finished.Mstep")
    }

    #### Check convergence of the algorithm and prepare for next iteration ####

    # Get the most promising theta value as given by the M-step
    mp.theta.next <- opt.res[which.max(opt.res[, "opt.EI"]), 1]

    # Check convergence
    converged <- EAM.converged(opt.theta.prev, evaluations, mp.theta.next,
                               iter.nbr, dir, hyperparams, verbose)

    # Increment iteration number
    iter.nbr <- iter.nbr + 1

    # Prepare for the next iteration: define vector of points to be checked in
    # the E-step of the next iteration.
    ptc.Estep <- c(opt.res[1], extra.points)
    opt.theta.prev <- theta.hash

  } # End main loop

  # Convergence of algorithm
  if (fps.fail_flag) {
    converge <- 2 # No initial feasible point was found
  } else if (iter.nbr <= hyperparams[["max.iter"]]) {
    converge <- 1 # Algorithm converged within max number of iterations
  } else {
    converge <- 0 # Algorithm didn't converge within max number of iterations
  }

  # Stop the chronometer
  if (time.run.duration) {
    chronometer$stop("algorithm finished")
  }

  # Return the results
  list(evaluations = evaluations, converge = converge,
       chronometer = chronometer)
}

#### Search strategies: main functions grid search ####

#' @title Set default hyperparameters for grid search algorithm
#'
#' @description This function returns a list with the (default) hyperparameters
#' used in the grid search algorithm
#'
#' @param options A list of user-specified values for (some of) the
#' hyperparameters. These hyperparameters could include:
#' \describe{
#'  \item{min.eval/max.eval:}{Minimum and maximum number of evaluations.}
#'  \item{next.gs.point:}{Function that determines the next point in the grid
#'  search sequence.}
#'  \item{step.size:}{Step size of the grid.}
#'  \item{bin.search.tol:}{Binary search tolerance.}
#'  \item{max.iter:}{Maximum number of iterations that the algorithm can run.}
#' }
#'
#' @returns List of hyperparameters for the gridsearch and binary search
#' algorithms.
#'
set.GS.hyperparameters <- function(options) {

  # Define the list of hyperparameters
  hyperparams <- list(

    # Minimum and maximum number of points to evaluate in initial feasible
    # search. min.eval must be at least 2 in order for the binary search
    # algorithm to work properly.
    min.eval = 10,
    max.eval = 100,

    # Type of grid search to be carried out
    next.gs.point = gs.binary,

    # Step size to be used in rudimentary grid search
    step.size = 1,

    # Convergence tolerance to be used in binary search grid search
    bin.search.tol = 1e-3,

    # Maximum number of iterations to run in the grid search algorithm
    max.iter = 100
  )

  # Precondition check
  if (identical(hyperparams[["next.gs.point"]], gs.binary) & (hyperparams[["min.eval"]] < 2)) {
    stop("When binary search is used, min.eval should be at least 2.")
  }

  # Overwrite default values with user specified values if needed
  hyperparams[names(options)] <- options

  # Return the result
  hyperparams
}

#' @title Grid search algorithm for finding the identified set
#'
#' @description This function implements the gridsearch and binary search
#' algorithms used to compute the roots of the violation curve and hence in
#' estimating the identified intervals.
#'
#' @param dir Search direction.
#' @param test.fun The test function to be inverted in order to obtain the
#' identified set. It should take a scalar parameter as argument (i.e. the
#' specified value of a component of the full parameter vector) and return a
#' list with named elements \code{list(theta, t.stat, crit.val)}, where
#' \code{theta} is the scalar value that was tested, \code{t.stat} is the value
#' of the test statistic and \code{crit.val} is the critical value to be used in
#' determining whether to reject or not reject.
#' @param hyperparams List of hyperparameters.
#' @param evaluations Matrix of already evaluated points, of which at least one
#' is feasible. When \code{evaluations = NULL} (default), the initial feasible
#' point search will be executed first.
#' @param time.run.duration Boolean value indicating whether to time each step
#' in the EAM algorithm. Requires \code{chronometer.R}. Default is
#' \code{time.run.duration = FALSE}.
#' @param verbose Boolean value indicating whether or not to print run time
#' updates to the console. Default is \code{verbose = FALSE}.
#' @param picturose Boolean value indicating whether or not to visualize the
#' identified set search. Default is \code{FALSE}.
#'
#' @returns List containing the evaluations of the test statistic and critical
#' values, convergence information, and run times.
#'
gridSearch <- function(dir, test.fun, hyperparams, evaluations = NULL,
                       time.run.duration = FALSE, verbose = 0,
                       picturose = FALSE) {

  # If the run times should be recorded, initialize a chronometer object
  if (time.run.duration) {
    source("chronometer.R")
    chronometer <- Chronometer$new()
    chronometer$start()
  } else {
    chronometer <- NULL
  }

  # Returns a set of points, of which at least one will be feasible
  if (!is.null(evaluations)) {
    fps.fail_flag <- FALSE
    if (all(evaluations[, "t.stat"] > evaluations[, "crit.val"])) {
      fps.fail_flag <- TRUE
    }
  } else {
    fps.out <- feasible_point_search(test.fun, hyperparams, verbose)
    evaluations <- fps.out[["evaluations"]]
    fps.fail_flag <- fps.out[["stop"]]
    if (time.run.duration) {
      chronometer$record.leg("finished.fps")
    }
  }

  # Initialize variables and function used in the main loop
  iter.nbr <- 0
  stop <- fps.fail_flag
  max.iter <- hyperparams[["max.iter"]]
  next.gs.point <- hyperparams[["next.gs.point"]]

  # Main loop
  if (verbose >= 3) {
    message("Entering main loop...")
  }
  while (!stop) {

    # Increment iteration number
    iter.nbr <- iter.nbr + 1

    # Get new r value and stopping criterion
    gnp.out <- next.gs.point(evaluations, dir, iter.nbr, hyperparams)
    r <- gnp.out[["theta.to.eval"]]
    stop <- gnp.out[["stop"]]

    # If a next test should be performed...
    if (!is.null(r)) {

      # Message to user. Update plot
      if (verbose >= 2) {
        message(sprintf("Iteration %s. Testing for r = %.3f", iter.nbr, r))
      }
      if (picturose) {
        plot_addpte(r)
      }

      # Perform the test
      res.Bei <- test.fun(r)

      # Store the results
      evaluations <- rbind(evaluations, res.Bei)

      # Record the run time of the iteration
      if (time.run.duration) {
        chronometer$record.leg(paste("Iter", iter.nbr))
      }

      # Update plot
      if (picturose) {
        plot_addpte.eval(evaluations[nrow(evaluations), , drop = FALSE])
      }

    } else if (!stop) {
      stop(paste("Next point is NULL while stop = FALSE. This should not",
                 "happen. Please contact the devs."))
    }
  } # End main loop

  # Convergence of algorithm
  if (fps.fail_flag) {
    converge <- 2 # No initial feasible point was found
  } else if (iter.nbr <= hyperparams[["max.iter"]]) {
    converge <- 1 # Algorithm converged within max number of iterations
  } else {
    converge <- 0 # Algorithm didn't converge within max number of iterations
  }

  # Stop the chronometer
  if (time.run.duration) {
    chronometer$stop("algorithm finished")
  }

  # Return the results
  list(evaluations = evaluations, chronometer = chronometer,
       converge = converge)
}





#### User interface: main estimation functions ####

#' @title Check argument consistency.
#'
#' @description This function checks whether the arguments supplied to the
#' main estimation function \code{pi.surv} are valid. When arguments are
#' invalid, the an exception is thrown.
#'
#' @inheritParams pi.surv
#'
check.args.pisurv <- function(data, idx.param.of.interest, idxs.c, t, par.space,
                              search.method, add.options) {

  #### Checks for the data ####

  # The provided data should be given as a data frame
  if (!is(data, "data.frame")) {
    stop("The provided data should be given as a data frame.")
  }

  # Check proper naming of columns
  colnames.unchecked <- colnames(data)
  if (!("Y" %in% colnames(data))) {
    stop(paste("The column containing the survival times should be named 'Y'",
               "(case sensitive)."))
  } else {
    colnames.unchecked <- setdiff(colnames.unchecked, "Y")
  }
  if (!("Delta" %in% colnames(data))) {
    stop(paste("The column containing the censoring indicator should be named",
               "Delta (case\n       sensitive)."))
  } else {
    colnames.unchecked <- setdiff(colnames.unchecked, "Delta")
  }
  if (!("X0" %in% colnames(data))) {
    stop(paste("The given data should contain an intercept column, named 'X0'",
               "(case sensitive)."))
  } else {
    colnames.unchecked <- setdiff(colnames.unchecked, "X0")
  }
  if (!("X1" %in% colnames(data))) {
    stop(paste("The given data should contain at least one covariate.",
               "covariates should be named\n       'X1', 'X2', ..."))
  }
  if (sum(grepl("[X][[:digit:]]+$", colnames(data))) !=
      max(as.numeric(gsub("X", "", colnames(data))[grepl("[X][[:digit:]]+$", colnames(data))])) + 1) {
    stop(paste("Invalid naming of the covariates detected. Covariates should",
               "be named 'X1', 'X2',\n       ... E.g. the case of two covariates",
               "named 'X1' and 'X3' is not allowed."))
  } else {
    colnames.unchecked <- setdiff(colnames.unchecked, colnames(data)[grepl("[X][[:digit:]]+$", colnames(data))])
  }
  if (length(colnames.unchecked) > 0) {
    stop(paste("Columns detected which do not have one of the",
               "required/allowed column names\n       (see documentation). Please",
               "remove these before excecution to ensure proper\n      ",
               "functioning of this function.",
               sprintf("\n\n       (Problematic columns: %s)",
                       paste(colnames.unchecked, collapse = ", "))))
  }

  # Check valid variable types
  if (any(apply(data, 2, class) != "numeric")) {
    stop("All variables should be given as numeric.")
  }

  # Check valid values for censoring indicator
  if (!all(data$Delta %in% c(0, 1))) {
    stop("Censoring indicator values can only be 0 or 1.")
  }

  # Check correct specification of idx.param.of.interest
  n.cov <- sum(grepl("X[[:digit:]]+$", colnames(data))) - 1
  n.param <- n.cov + 1
  if (is(idx.param.of.interest, "character")) {
    if (idx.param.of.interest != "all") {
      stop("Invalid specification for idx.param.of.interest.")
    }
  } else {
    if (!(class(idx.param.of.interest)) %in% c("numeric", "integer")) {
      stop("When idx.param.of.interest is not 'all', it should be an integer")
    }
    if (idx.param.of.interest != round(idx.param.of.interest)) {
      stop("When idx.param.of.interest is not 'all', it should be an integer")
    }
    if (idx.param.of.interest > n.param | idx.param.of.interest < 1) {
      stop("Invalid index of parameter of interest")
    }
  }

  # Check correct specification of idxs.c
  if (!all(class(idxs.c) %in% c("numeric", "integer"))) {
    stop("idxs.c should be a vector of numerics.")
  }
  if (!all(idxs.c == round(idxs.c))) {
    stop("idxs.c should be integer valued.")
  }
  if (any(idxs.c > n.cov) | any (idxs.c < 1)) {
    stop("Invalid indices detected in idxs.c")
  }

  # Check valid specification of time point of interest
  if (!is(t, "numeric")) {
    stop("Time point of interest t should be numeric.")
  }

  # Check correct specification of parameter space.
  if (!("matrix" %in% class(par.space))) {
    stop("Class of par.space should be 'matrix'.")
  }
  if (any(apply(par.space, 1:2, class) != "numeric")) {
    stop(paste("Elements of the matrix containing the bounds on the parameter",
               "space should be numeric."))
  }
  if (ncol(par.space) != 2) {
    stop("par.space must have 2 columns.")
  }
  if (any(apply(par.space, 1, function(row) {row[1] >= row[2]}))) {
    stop("Invalid bounds on parameter space detected.")
  }
  if (nrow(par.space) != n.cov + 1) {
    stop("Invalid number of rows in parameter space given the covariates",
         "provided in the data set.")
  }

  # Check correct specification of search.method.
  if (!is(search.method, "character")) {
    stop("Invalid specification of search method.")
  }
  if (!(search.method %in% c("EAM", "GS"))) {
    stop("Invalid specification of search.method")
  }
}

#' @title Partially identify the coefficients in the model
#' \eqn{\Lambda(x^\top \beta(t))} for the given data set. This methodology
#' implements the one described in Willems et al. (2024+).
#'
#' @description This function estimates bounds on the coefficients the single-
#' index model \eqn{\Lambda(x^\top \beta(t))} for the conditional CDF of the
#' event time.
#'
#' @param data Data frame containing the data on which to fit the model. The
#' columns should be named as follows: 'Y' = observed timed, 'Delta' = censoring
#' indicators, 'X0' = intercept column, 'X1' - 'Xp' = covariates.
#' @param idx.param.of.interest Index of element in the covariate vector for
#' which the identified interval should be estimated. It can also be specified
#' as \code{idx.param.of.interest = "all"}, in which case identified intervals
#' will be computed for all elements in the parameter vector.
#' @param idxs.c Vector of indices of the continuous covariates. Suppose the
#' given data contains 5 covariates, of which 'X2' and 'X5' are continuous, this
#' argument should be specified as \code{idxs.c = c(2, 5)}.
#' @param t Time point for which to estimate the identified set of
#' \eqn{\beta(t)}.
#' @param par.space Matrix containing bounds on the space of the parameters. The
#' first column corresponds to lower bounds, the second to upper bounds. The i'th
#' row corresponds to the bounds on the i'th element in the parameter vector.
#' @param search.method The search method to be used to find the identified
#' interval. Default is \code{search.method = "GS"}.
#' @param add.options List of additional options to be specified to the method.
#' These options can range from 'standard' hyperparameters such as the
#' confidence level of the test and number of instrumental functions to be used,
#' to technical hyperparameters regarding the search method and test
#' implementation. For the latter, we refer to the documentations of
#' \code{set.hyperparameters}, \code{set.EAM.hyperparameters} and
#' \code{set.GS.hyperparameters}. We recommend to use the default parameters,
#' unless you really know what you are doing.
#' @param picturose Picturosity flag. If \code{TRUE}, a plot illustrating the
#' workings of the algorithm will updated during runtime. Default is
#' \code{picturose = FALSE}.
#' @param parallel Flag for whether or not parallel computing should be used.
#' Default is \code{parallel = FALSE}. When \code{parallel = TRUE}, this
#' implementation will use \code{min(detectCores() - 1, 10)} cores to construct
#' the parallel back-end.
#'
#' @returns Matrix containing the identified intervals of the specified
#' coefficients, as well as corresponding convergence information of the
#' estimation algorithm.
#'
#' @export
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#'
#' @examples
#' \donttest{
#'   #     - Link function: AFT link function (default setting)
#'   #     - Number of IF: 5 IF per continuous covariate (default setting)
#'   #     - Search method: Binary search
#'   #     - Type of IF: Cubic spline functions for continuous covariate, indicator
#'   #       function for discrete covariate (default setting).
#'
#'   # Load the 'depCensoring' package
#'   library("depCensoring")
#'
#'   # Load 'survival' package in R.
#'   library("survival")
#'
#'   # Load and preprocess data
#'   data <- survival::lung
#'   data[, "intercept"] <- rep(1, nrow(data))
#'   data[, "status"] <- data[, "status"] - 1
#'   data <- data[, c("time", "status", "intercept", "age", "sex")]
#'   colnames(data) <- c("Y", "Delta", "X0", "X1", "X2")
#'
#'   # Settings for main estimation function
#'   idx.param.of.interest <- 1 # Interest in effect of age
#'   idxs.c <- 1                # X1 (age) is continuous
#'   t <- 200                   # Model imposed at t = 200
#'   search.method <- "GS"      # Use binary search
#'   par.space <- matrix(rep(c(-10, 10), 3), nrow = 3, byrow = TRUE)
#'
#'   # Estimate the identified intervals
#'   pi.surv(data, idx.param.of.interest, idxs.c, t, par.space, search.method)
#' }
#'
#'
#' @references Willems, I., Beyhum, J. and Van Keilegom, I. (2024+). Partial
#' identification for a class of survival models under dependent censoring.
#' (In preparation).
#'
pi.surv <- function(data, idx.param.of.interest, idxs.c, t, par.space,
                    search.method = "GS", add.options = list(),
                    picturose = FALSE, parallel = FALSE) {

  #### Consistency checks ####

  check.args.pisurv(data, idx.param.of.interest, idxs.c, t, par.space,
                    search.method, add.options)

  # If required, set-up parallel back-end.
  if (parallel) {
    n.cores <- min(detectCores() - 1, 10)
    clust <- makeCluster(n.cores)
    registerDoParallel(clust)
    add.options$clust <- clust
  }

  #### Set the hyperparameters of the algorithm ####

  # Number of covariates
  n.cov <- sum(grepl("X[1-9][[:digit:]]*", colnames(data)))

  # Transform argument 'idx.param.of.interest' to list of unit vectors with a
  # 1 placed at the respective index.
  c.to.check <- list()
  if (is(idx.param.of.interest, "character")) {
    for (i in 1:(n.cov + 1)) {
      c.vec <- rep(0, n.cov + 1)
      c.vec[i] <- 1
      c.to.check <- c(c.to.check, list(c.vec))
    }
  } else {
    c.vec <- rep(0, n.cov + 1)
    c.vec[idx.param.of.interest] <- 1
    c.to.check <- list(c.vec)
  }

  # Pre-set algorithm hyperparameters
  options <- list(n.if.per.cov = 5,
                  K.bar = 3,
                  B = 600,
                  next.gs.point = gs.binary,
                  alpha = 0.95,
                  link.function = "AFT_ll",
                  inst.func.family = "cd",
                  degree = 3,
                  G.c = G.spline,
                  idxs.c = idxs.c)

  # Overwrite pre-sets with user-specified options
  options[names(add.options)] <- add.options

  #### Find the identified interval(s) ####

  # Initialize object that will store the results
  results <- matrix(nrow = length(c.to.check), ncol = 4)
  rownames(results) <- paste0("X", matrix(unlist(c.to.check), nrow = length(c.to.check)) %*% 0:n.cov)
  colnames(results) <- c("lower", "upper", "conv.l", "conv.u")

  # Run the algorithm for the selected elements of the parameter vector
  for (c in c.to.check) {

    # Get the identified set of the covariate of interest
    fis.out <- find.identified.set(c, t, par.space, data, search.method, options,
                                   verbose = 3, picturose = picturose,
                                   parallel = parallel)

    # Store the results
    row.name <- paste0("X", which(c == 1) - 1)
    results[row.name, c("lower", "upper")] <- fis.out$ident.set
    results[row.name, c("conv.l", "conv.u")] <- c(fis.out$converge1, fis.out$converge2)
  }

  # Remove parallel back-end
  if (parallel) {
    stopCluster(clust)
  }

  # Return the results
  results
}

#' @title Subset a given data set based on the provided criteria.
#'
#' @description This function subsets a data set based on the given criteria.
#'
#' @param data Data frame (see documentation of 'pi.surv.R')
#' @param criteria A data frame wherein each row contains the variable for which
#' the criterion should apply, a comparing operator, and the reference value.
#' @param keep.covs Names of covariates to keep. Default \code{keep.covs = NULL}.
#' @param max.row Maximum rows in the resulting data set. Default is
#' \code{max.row = NULL}
#'
#' @noRd
#'
subset.data <- function(data, criteria, keep.covs = NULL, max.row = NULL) {
  data.sub <- data
  for (row.idx in 1:nrow(criteria)) {

    # Get variable name or index for which to apply the criterion
    var.name <- criteria[row.idx, 1]

    # Apply criterion over data set
    expr <- sprintf("data.sub[which(data.sub[, var.name] %s criteria[row.idx, 3]),]",
                    criteria[row.idx, 2])
    data.sub <- eval(parse(text = expr))
  }

  # Retain only the specified covariates
  name.dict <- NULL
  if (!is.null(keep.covs)) {
    data.sub <- data.sub[, c("Y", "Delta", "X0", keep.covs)]
    name.dict <- data.frame(
      old = c("Y", "Delta", "X0", keep.covs),
      new = c("Y", "Delta", paste0("X", 0:length(keep.covs)))
    )
    colnames(data.sub) <- name.dict[match(colnames(data.sub), name.dict[, 1]), 2]
  }

  # Retain only the first 'max.row' rows of the data frame
  if (!is.null(max.row)) {
    data.sub <- data.sub[1:min(nrow(data.sub), max.row), ]
  }

  # Return subsetted data and additional information
  list(data.sub, name.dict)
}

#' @title Combine bounds based on majority vote.
#'
#' @description This function combines a list of individual identified intervals
#' to a single one based on majority vote. Note that the intersection of all
#' intervals can be viewed as a majority vote as well, so that it is included as
#' a special case.
#'
#' @param results.list List object containing the individual identified
#' intervals.
#' @param threshold Threshold proportion of identified intervals a given value
#' should be contained in in order for it to be included in the combined
#' identified interval. For intersection bounds, set this value to \code{1}.
#'
#' @returns The combined identified interval.
#'
cbMV <- function(results.list, threshold) {

  # Obtain some information about the results
  cov.names <- rownames(results.list[[1]])

  # Obtain vector of lower and upper bounds
  lbs <- matrix(ncol = length(cov.names), nrow = length(results.list))
  ubs <- matrix(ncol = length(cov.names), nrow = length(results.list))
  missspec <- FALSE
  for (results.idx in 1:length(results.list)) {

    # Get results of this iteration
    results <- results.list[[results.idx]]

    # Extract the bounds on the parameters
    for (cov.name.idx in 1:length(cov.names)) {
      cov.name <- cov.names[cov.name.idx]
      lbs[results.idx, cov.name.idx] <- results[cov.name, "lower"]
      ubs[results.idx, cov.name.idx] <- results[cov.name, "upper"]
    }

    # Check if the model was determined to be misspecified.
    if (any(results[, c("conv.l", "conv.u")] == 2)) {
      missspec <- TRUE
    }
  }
  colnames(lbs) <- colnames(ubs) <- cov.names

  # Initialize object that will store the combined bounds
  bounds.MV <- matrix(rep(c(-Inf, Inf, 1, 1), length(cov.names)),
                      nrow = length(cov.names),
                      byrow = TRUE)
  colnames(bounds.MV) <- colnames(results.list[[1]])
  rownames(bounds.MV) <- cov.names

  # If model was misspecified at any of the tested points, also the combined
  # model is misspecified.
  if (missspec) {
    bounds.MV[, 3:4] <- 2
    return(bounds.MV)
  }

  # Else, determine the combined bounds for each parameter of interest.
  for (cov.name in cov.names) {

    # Get bounds of each part
    ths <- sort(unique(c(lbs[, cov.name], ubs[, cov.name])))
    parts <- matrix(rep(ths, c(1, rep(2, length(ths) - 2), 1)), ncol = 2, byrow = TRUE)

    # For each part, obtain the number of votes
    parts <- cbind(parts, rep(0, nrow(parts)))
    for (part.idx in 1:nrow(parts)) {
      parts[part.idx, 3] <- sum((lbs[, cov.name] < mean(parts[part.idx, 1:2])) &
                                  (ubs[, cov.name] > mean(parts[part.idx, 1:2])))
    }

    # Subset to parts getting majority vote
    MV.parts <- parts[parts[, 3]/length(lbs[, cov.name]) >= threshold, 1:2, drop = FALSE]

    # Recombine neighbouring bounds
    row.idx <- 1
    while (row.idx < nrow(MV.parts)) {
      if (MV.parts[row.idx, 2] == MV.parts[row.idx + 1, 1]) {
        MV.parts[row.idx, 2] <- MV.parts[row.idx + 1, 2]
        MV.parts <- MV.parts[setdiff(1:nrow(MV.parts), row.idx + 1), , drop = FALSE]
      } else {
        row.idx <- row.idx + 1
      }
    }

    # If the result is not an interval, take lower and upper bounds, but warn
    # the user.
    if (nrow(MV.parts) > 1) {
      MV.parts <- matrix(c(min(MV.parts), max(MV.parts)), nrow = 1)
      warning("Non-connected identified set. Returning outer identified interval.")
    }

    # Store the result
    bounds.MV[cov.name, ] <- cbind(MV.parts, 1, 1)
  }

  # Return the result
  bounds.MV
}

#' @title Get the file name of the results from Pancreas data applications
#'
#' @description This function obtains the names of the files containing the
#' results of the Pancreas data application.
#'
#' @param args The vector of arguments supplied to the data application wrapper
#' function.
#'
#' @noRd
#'
get.file.name.pancreas <- function(args, idx.param.of.interest, master.dir) {

  # Create directory if necessary
  check_create.dir(master.dir)

  # Get name of variable
  if (idx.param.of.interest == 1) {
    var.name <- "intercept"
  } else if (idx.param.of.interest == 2) {
    var.name <- "age"
  } else if (idx.param.of.interest == 3) {
    var.name <- "tumorsize"
  }

  # Create file name
  file.name <- sprintf("crit-%d__nifpc-%d_var-%s.csv", args["criteria.idx"],
                       args["n.if.per.cov"], var.name)

  # Return path name
  paste0(c(master.dir, file.name), collapse = "/")
}

#### User interface: plotting functions (picturose) ####

#' @title Clear plotting window
#'
#' @description This function clears the plotting window
#'
#' @importFrom grDevices dev.off
#'
clear.plt.wdw <- function() {
  tryCatch(invisible(dev.off()), error = function(e) {invisible(e)})
}

#' @title Draw base plot
#'
#' @description This functon draws the base plot, used when
#' \code{picturose = TRUE}.
#'
#' @param c Projection vector
#' @param hp List of hyperparameters
#'
plot_base <- function(c, hp) {

  # Clear plotting window
  clear.plt.wdw()

  # Make base plot
  c.idx <- which(c == 1) - 1
  plot(x = c(hp$theta.lb, hp$theta.ub), y = c(0, 0), type = 'l',
       xlab = bquote("parameter"~"space"~"for"~theta[.(c.idx)]),
       ylab = "",
       yaxt = 'n')
}

#' @title Draw points to be evaluated
#'
#' @description This function draws the points to be evaluated.
#'
#' @param pte Vector of points to be evaluated.
#' @param col Color of the points.
#'
#' @importFrom graphics points
#'
plot_addpte <- function(pte, col = "orange") {
  points(x = pte, y = rep(0, length(pte)), pch = 16, col = col)
  points(x = pte, y = rep(0, length(pte)), pch = 1, col = "black")
}

#' @title Draw evaluated points.
#'
#' @description This function draws evaluated points. Feasible points are
#' indicated in green, red points correspond to infeasible points.
#'
#' @param evaluations Matrix of evaluations to be drawn.
#'
#' @importFrom graphics points
#'
plot_addpte.eval <- function(evaluations) {
  feas <- evaluations[, 2] <= evaluations[, 3]
  col <- ifelse(feas, "green", "red")
  points(x = evaluations[, 1], y = rep(0, nrow(evaluations)), col = col, pch = 16)
  points(x = evaluations[, 1], y = rep(0, nrow(evaluations)), pch = 1, col = "black")
}






#### Simulation functions ####

#' @title Define the hyperparameters used for finding the identified interval
#'
#' @description This function defines all the necessary hyperparameters used to
#' run the methodology.
#'
#' @param data Data frame.
#' @param par.space Bounds on the parameter space.
#' @param c Projection vector.
#' @param search.method Search method to use (\code{"EAM"} or \code{"GS"})
#' @param options List of user specified hyperparameters that will substitute
#' the corresponding default values. This list can contain the entries:
#' \describe{
#'  \item{cov.ranges:}{known bounds on each of the covariates in the data set.}
#'  \item{norm.func.name:}{Name of the normalization function to be used. Can
#'  be either "normalize.covariates1" or "normalize.covariates2" (default).
#'  The former is a simple elementwise rescaling. The latter uses the PCA
#'  approach as discussed in Willems et al. (2024+).}
#'  \item{inst.func.family:}{Family of instrumental functions to be used for
#'  all covariates. Options are "box", "spline" and "cd". The former two are
#'  only applicable for continuous covariates. The latter can also handle
#'  discrete covariates. Default is "cd".}
#'  \item{G.c:}{The class of instrumental functions used for the continuous
#'  covariates in the model, in case "cd" is selected as
#'  \code{inst.func.family:}. Options are "box" and "spline". Default is
#'  "spline".}
#'  \item{degree:}{The degree of the B-spline functions, should they be used as
#'  instrumental functions for the continuous covariates. Default is 3.}
#'  \item{link.function:}{Name of the link function to be used. Options are
#'  "AFT_ll" for the AFT model with log-logistic baseline, or "Cox_wb" for the
#'  Cox PH model (originally with Weibull baseline, but now for a general)
#'  baseline hazard).}
#'  \item{K.bar:}{Number of refinement steps when obtaining the critical value.
#'  See Bei (2024).}
#'  \item{B:}{Number of bootstrap samples to be used when obtaining the
#'  bootstrap distribution of the test statistic.}
#'  \item{ignore.empty.IF:}{Boolean value indicating whether instrumental
#'  functions with empty support should be ignored (cf. Willems et al., 2024).
#'  Default is FALSE. The feature \code{ignore.empty.IF = TRUE} is experimental,
#'  so there might exist edge cases for which the implementation will fail to
#'  run.}
#' }
#' Other (hidden) options can also be overwritten, though we highly discourage
#' this. If necessary, you can consult the source code of this functions to
#' find the names of the desired parameters and add their name alongside their
#' desired value as an entry in \code{options} (e.g.
#' \code{options$min.var <- 1e-4}. Again, not recommended!).
#'
#' @returns The list of hyperparameters.
#'
set.hyperparameters <- function(data, par.space, c, search.method, options) {

  #### General tuning hyperparameters ####

  # Number of covariates (excluding the intercept)
  n.cov <- sum(grepl("X[[:digit:]]+", colnames(data))) - 1
  n <- nrow(data)

  # Type of covariates. DGP option is used when running the simulations of
  # Willems et al. (2024+).
  cov.idxs <- 1:n.cov
  if (is.null(options[["idxs.c"]])) {
    if (options[["DGP"]] <= 20) {
      idxs.c <- cov.idxs
    } else if (20 < options[["DGP"]] & options[["DGP"]] <= 40) {
      idxs.c <- cov.idxs[1:ceiling(n.cov/2)]
    } else if (40 < options[["DGP"]] & options[["DGP"]] <= 60) {
      idxs.c <- integer(0)
    }
  } else {
    idxs.c <- options[["idxs.c"]]
  }
  if (!all(idxs.c %in% 1:n.cov)) {
    stop("Invalid value for 'idxs.c'")
  }

  # Tuning parameters
  kappa.n <- sqrt(log(n))
  lambda.n <- sqrt(n) * kappa.n^2
  epsilon.n <- sqrt(log(kappa.n^2)/n)
  delta.n <- min(1/n, 10^(-4))

  # Artificial minimum variance (ensure non-nullity of variance).
  min.var <- 1e-6

  # Instrumental function hyperparameters
  n.if.per.cov <- options[["n.if.per.cov"]]
  if (is.null(n.if.per.cov)) {
    stop("Number of instumental functions per covariate should be specified")
  }
  inst.func.family <- options[["inst.func.family"]]
  if (is.null(inst.func.family)) {
    inst.func.family <- "cd"
  }
  ignore.empty.IF <- FALSE
  if ("ignore.empty.IF" %in% names(options)) {
    if (is(options[["ignore.empty.IF"]], "logical")) {
      ignore.empty.IF <- options[["ignore.empty.IF"]]
    } else {
      stop("Specified value for 'ignore.empty.IF' should be TRUE/FALSE.")
    }
  }
  cov.ranges <- options[["cov.ranges"]]

  # Normalization function (used in instrumental function)
  norm.func.name <- options[["norm.func.name"]]
  if (is.null(norm.func.name)) {
    norm.func.name <- "normalize.covariates2"
  }
  if (norm.func.name == "normalize.covariates1") {
    norm.func <- normalize.covariates
  } else if (norm.func.name == "normalize.covariates2") {
    norm.func <- normalize.covariates2
  }

  # Precompute normalized covariates. When inst.func.family == cd, this should
  # only be done on the continuous covariates.
  if (inst.func.family == "cd") {
    cols.to.include <- which(!(colnames(data) %in% paste0("X", setdiff(1:n.cov, idxs.c))))
    data.c <- data[, cols.to.include]
    cov.ranges.c <- cov.ranges[, cols.to.include]
    norm.cov.out <- norm.func(data = data.c, cov.ranges = cov.ranges.c,
                              idxs.c = "all")
  } else {
    norm.cov.out <- norm.func(data = data, cov.ranges = cov.ranges)
  }

  # Number of instrumental functions (computed differently depending on whether
  # or not G.cd or G.cd.mc is selected).
  if (inst.func.family == "cd") {

    # Number of instrumental functions pertaining to continuous covariates
    names.cov.d <- setdiff(paste("X", setdiff(1:n.cov, idxs.c), sep = ""), "X")
    n.inst.func.c <- n.if.per.cov^(n.cov - length(names.cov.d))

    # Number of instrumental functions pertaining to discrete covariates, taking
    # into account that some combinations of discrete covariate levels might
    # be empty.
    covariates.d <- data[, names.cov.d, drop = FALSE]
    n.inst.func.d <- max(nrow(unique(covariates.d)), 1)

    # Total number of covariates
    n.inst.func <- n.inst.func.c * n.inst.func.d

  } else if (inst.func.family == "cd.manycov") {

    # Precompute some necessary information
    cov.names <- colnames(data)[grep("X[1-9][[:digit:]]*$", colnames(data))]
    info.manycov <- data.frame(cov.pair = character(), n.if = numeric())
    for (cov.name.idx1 in 1:length(cov.names)) {
      for (cov.name.idx2 in 2:length(cov.names)) {
        if (cov.name.idx2 > cov.name.idx1) {

          # Name of covariates in the pair
          cov.name1 <- cov.names[cov.name.idx1]
          cov.name2 <- cov.names[cov.name.idx2]

          # Number of instrumental functions for each
          n.if1 <- ifelse(cov.name.idx1 %in% idxs.c, n.if.per.cov, length(unique(data[, cov.name1])))
          n.if2 <- ifelse(cov.name.idx2 %in% idxs.c, n.if.per.cov, length(unique(data[, cov.name2])))

          # Total number of instrumental functions
          n.if <- n.if1 * n.if2

          # Add to information data frame
          row <- list(cov.pair = sprintf("%s, %s", cov.name1, cov.name2),
                      n.if = n.if)
          info.manycov <- rbind(info.manycov, row)
        }
      }
    }

    # Add supplementary rows and columns
    info.manycov <- cbind(idx = 1:nrow(info.manycov),
                          info.manycov,
                          cumsum = cumsum(info.manycov$n.if))
    info.manycov <- rbind(list(idx = 0, cov.pair = "", n.if = 0, cumsum = 0),
                          info.manycov)

    # Get number of instrumental functions
    n.inst.func <- max(info.manycov$cumsum)

  } else {
    n.inst.func <- n.if.per.cov^n.cov
  }

  # If "G.cd" is selected, obtain all levels of the 'combined' discrete
  # covariates present in the data.
  discrete.covariate.levels <- NULL
  if (inst.func.family == "cd") {
    names.cov.d <- setdiff(paste("X", setdiff(1:n.cov, idxs.c), sep = ""), "X")
    discrete.covariate.levels <- unique(data[, names.cov.d, drop = FALSE])
  }

  # Instrumental function
  if (n.cov > 3 & inst.func.family != "cd.manycov") {
    warning(paste("It is recommended to specify 'inst.func.family = 'cd.manycov'",
                  "when the number of covariates exceeds 3."))
  }
  if (inst.func.family == "box") {
    degree <- 0
    G <- function(x, g.idx) {G.box(x, g.idx, data, n.if.per.cov,
                                   cov.ranges = cov.ranges,
                                   norm.func = norm.func,
                                   norm.cov.out = norm.cov.out)}

  } else if (inst.func.family == "spline") {
    degree <- options[["degree"]]
    if (is.null(degree)) {degree <- 3}
    G <- function(x, g.idx) {G.spline(x, g.idx, data, n.if.per.cov,
                                      degree = degree, cov.ranges = cov.ranges,
                                      norm.func = norm.func,
                                      norm.cov.out = norm.cov.out)}

  } else if (inst.func.family == "cd") {
    G.c <- options[["G.c"]]
    if (is.null(G.c)) {stop("G.c should be specified when G.cd is used")}
    degree <- options[["degree"]]
    if (is.null(degree)) {degree <- 3}
    G <- function(x, g.idx) {G.cd(x = x, g.idx = g.idx, data = data,
                                  n.if.per.cov = n.if.per.cov, idxs.c = idxs.c,
                                  G.c = G.c, norm.func = norm.func,
                                  discrete.covariate.levels = discrete.covariate.levels,
                                  cov.ranges = cov.ranges,
                                  norm.cov.out = norm.cov.out, degree = degree)}

  } else if (inst.func.family == "cd.manycov") {
    G.c <- options[["G.c"]]
    if (is.null(G.c)) {stop("G.c should be specified when G.cd is used")}
    degree <- options[["degree"]]
    if (is.null(degree)) {degree <- 3}
    G <- function(x, g.idx) {
      G.cd.mc(x = x, g.idx = g.idx, data = data, n.if.per.cov = n.if.per.cov,
              idxs.c = idxs.c, G.c = G.c, norm.func = norm.func,
              info.manycov = info.manycov, cov.ranges = cov.ranges,
              degree = degree)
    }

  } else {
    stop("Unknown instrumental function family specified.")
  }

  # Link function to be used
  link.function <- "AFT_ll"
  if ("link.function" %in% names(options)) {
    link.function <- options[["link.function"]]
  }
  if (link.function == "AFT_ll") {
    Lambda <- Lambda_AFT_ll
    dLambda <- dLambda_AFT_ll
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (link.function == "Cox_wb") {
    Lambda <- Lambda_Cox_wb
    dLambda <- dLambda_Cox_wb
    inv.Lambda <- Lambda_inverse_Cox_wb
  }

  # Number of initial values in step 4 of the algorithm of Bei (2024)
  K.bar <- options[["K.bar"]]
  if (is.null(K.bar)) {
    stop(paste0("Number of initial values in step 4 of the algorithm of Bei ",
                "(2024) should be specified."))
  }

  # Number of bootstrap samples
  B <- options[["B"]]
  if (is.null(B)) {
    stop("Number of bootstrap samples should be specified")
  }

  hp <- list(

    # Link function
    "Lambda" = Lambda,
    "dLambda" = dLambda,
    "inv.Lambda" = inv.Lambda,

    # Instrumental functions
    "n.if.per.cov" = n.if.per.cov,
    "G" = G,
    "n.inst.func" = n.inst.func,
    "discrete.covariate.levels" = discrete.covariate.levels,
    "norm.func.name" = norm.func.name,

    # Tuning parameters
    "kappa.n" = kappa.n,
    "lambda.n" = lambda.n,
    "epsilon.n" = epsilon.n,
    "delta.n" = delta.n,
    "min.var" = min.var,

    # Bounds of the parameter space for theta
    "theta.lb" = par.space[which(c == 1), 1],
    "theta.ub" = par.space[which(c == 1), 2],

    # Full parameter space
    "par.space" = par.space,

    # Hyperparameters for computing the test statistic and critical value
    "K.bar" = K.bar,    # Number of initial values in step 4 of algorithm
    # described in Bei (2024).
    "B" = B             # Number of bootstrap samples
  )

  # Overwrite default values with user specified values if needed
  hp[names(options)] <- options

  #### Search strategy-specific hyperparameters ####

  if (search.method == "EAM") {
    ss.hp <- set.EAM.hyperparameters(options)
  } else if (search.method == "GS") {
    ss.hp <- set.GS.hyperparameters(options)
  }

  #### Return all hyperparamters ####

  hp[names(ss.hp)] <- ss.hp
  return(hp)
}

#' @title Estimate the identified set of the parameters in the model.
#'
#' @description This function estimates the identified set of the parameters in
#' the model. It does so by running the selected search algorithm twice (once
#' for minimization of the nonlinear program and once for maximization). Not to
#' be confused with the low level function 'get.identified.set.R'.
#'
#' @param c Projection vector
#' @param t Time point of interest.
#' @param par.space Bounds on the parameter space.
#' @param data Data frame.
#' @param search.method String value indicating the search method to be used
#' for finding the identified set. Can be 'EAM' or 'GS' (grid search).
#' @param options List of user specified hyperparameters that will substitute
#' the corresponding default values.
#' @param verbose Verbosity parameter. Higher values indicate larger degrees of
#' verbosity. Default is \code{verbose = 0} (no verbosity).
#' @param time.run.duration Boolean value indicating whether run durations
#' should be timed.
#' @param picturose Picturosity flag. If \code{TRUE}, a plot illustrating the
#' workings of the algorithm will updated during runtime. Default is
#' \code{picturose = FALSE}.
#' @param parallel Flag for whether or not parallel computing should be used.
#' Default is \code{parallel = FALSE}.
#'
#' @returns List containing the estimated identified set, convergence
#' information and run times.
#'
#' @noRd
#'
find.identified.set <- function(c, t, par.space, data, search.method, options,
                                verbose = 0, time.run.duration = FALSE,
                                picturose = FALSE, parallel = FALSE) {

  #### Precondition checks ####

  if (verbose >= 2) {
    message("  Performing precondition checks...")
  }

  # Check if the projection vector is valid
  if (!all(c %in% c(0, 1)) | sum(c^2) != 1) {
    stop("Invalid argument for c")
  }

  # Check if the bounds of the parameter space are valid
  if (!is.matrix(par.space)) {
    stop("par.space should be a matrix")
  } else {
    if (ncol(par.space) != 2) {
      stop("par.space should contain 2 columns")
    }
  }

  # Check if the data is specified as required
  if (!("X0" %in% colnames(data))) {
    stop("Data frame should contain intercept column named 'X0'.")
  }
  if (!any(grepl("X[1-9]+", colnames(data)))) {
    warning("No covariates detected in the data")
  }

  #### Initialize hyperparameters ####

  if (verbose >= 2) {
    message("  Initializing hyperparameters...")
  }

  # Hyperparameters
  hp <- set.hyperparameters(data, par.space, c, search.method, options)

  # Set search algorithm
  if (search.method == "EAM") {
    search.algorithm <- EAM
  } else if (search.method == "GS") {
    search.algorithm <- gridSearch
  }

  # Precompute instrumental function evaluations
  inst.func.evals <- t(get.instrumental.function.evals(data, hp))

  # Find possible instrumental functions with empty support. If indicated,
  # remove these from the analysis.
  idx.empty.IF <- which(rowSums(inst.func.evals) == 0)
  if (!is.null(hp$ignore.empty.IF)) {
    if (hp$ignore.empty.IF & length(idx.empty.IF) > 0) {

      # Remove empty IF
      inst.func.evals <- inst.func.evals[setdiff(1:nrow(inst.func.evals), idx.empty.IF), ]

      # Update the hyperparameter information
      hp$n.inst.func <- nrow(inst.func.evals)

      # Warn the user
      warning(sprintf("%d empty instrumental functions removed!", length(idx.empty.IF)))
    }
  }

  # Test whether each instrumental function contains at least one observation.
  # If not, return with appropriate failure flag
  if (any(rowSums(inst.func.evals) == 0)) {
    return(list(ident.set = c(-Inf, Inf),
                converge1 = 3,
                converge2 = 3,
                chronometer1 = Chronometer$new(),
                chronometer2 = Chronometer$new()))
  }

  # Set test function
  test.fun <- function(theta) {
    test.point_Bei(theta, c, t, par.space, data, hp, verbose = FALSE,
                   inst.func.evals = inst.func.evals, alpha = options$alpha)
  }

  #### Find identified set ####

  # If required, update user and initialize plot
  if (verbose >= 2) {
    message("  Starting search for initial feasible points...")
  }
  if (picturose) {
    plot_base(c, hp)
  }

  # Pre-search for feasible points
  fps.out <- feasible_point_search(test.fun, hp, verbose, picturose = picturose,
                                   parallel = parallel)
  evaluations <- fps.out[["evaluations"]]

  # Re-set test function, this time using parallel computing, if necessary.
  test.fun <- function(theta) {
    test.point_Bei(theta, c, t, par.space, data, hp, verbose = FALSE,
                   inst.func.evals = inst.func.evals, alpha = options$alpha,
                   parallel = parallel)
  }

  # Run search algorithm in dir = 1
  if (verbose >= 2) {
    message("  Starting search in dir = 1")
  }
  dir <- 1
  sa.out <- search.algorithm(dir, test.fun, hp, evaluations, time.run.duration,
                             verbose, picturose = picturose)
  evaluations1 <- sa.out[["evaluations"]]
  converge1 <- sa.out[["converge"]]
  chronometer1 <- sa.out[["chronometer"]]

  # Run search algorithm in dir = -1
  if (verbose >= 2) {
    message("  Starting search in dir = -1")
  }
  dir <- -1
  sa.out <- search.algorithm(dir, test.fun, hp, evaluations, time.run.duration,
                             verbose, picturose = picturose)
  evaluations2 <- sa.out[["evaluations"]]
  converge2 <- sa.out[["converge"]]
  chronometer2 <- sa.out[["chronometer"]]

  # Combine results and find identified interval
  evaluations <- rbind(evaluations1, evaluations2)
  ident.set <- get.identified.set(evaluations)

  #### Return the results ####

  list(ident.set = ident.set, converge1 = converge1, converge2 = converge2,
       chronometer1 = chronometer1, chronometer2 = chronometer2)
}

#' @title Obtain the file name of a simulation
#'
#' @description This function returns the filename, prefixed with relevant
#' folder paths, for the results of a simulation iteration. (This function is
#' used in simulate1D.R)
#'
#' @param comb Either a matrix or list containing the information of the
#' simulation design.
#' @param sim.iter.nbr Number of the simulation ran (value in 1:n.sim)
#' @param seed The random seed for which the simulation is ran.
#' @param master.dir Master directory.
#' @param shortened Boolean value indicating whether a shortened version of the
#' file name should be returned.
#' @param defaults Default values of parameters that will be omitted from file
#' name if \code{shortened = TRUE}.
#'
#' @note This function creates the relevant directories if necessary by calling
#' the function 'lowLevelFunctions::check_create.dir.R'.
#'
#' @noRd
#'
get.file.name <- function(comb, sim.iter.nbr, seed, master.dir,
                          shortened = FALSE, defaults = NULL) {

  if (is(comb, "matrix")) {

    # Shorten inst.func.family name, if applicable.
    if (comb[, "inst.func.family"] == "cd.manycov") {
      comb[, "inst.func.family"] <- "cdmc"
    }

    # Directory: (link function,) search method, inst.func.family, K.bar
    sub.dir <- sprintf("Search-%s__IF-%s__Kbar-%s__Alpha-%s", comb[, "search.method"],
                       comb[, "inst.func.family"], comb[, "K.bar"], comb[, "alpha"])
    if ("link.function" %in% colnames(comb)) {
      sub.dir <- sprintf("LF-%s__%s", comb[, "link.function"], sub.dir)
    }
    if ("parametric" %in% colnames(comb)) { # Option used in misspecification simulations
      sub.dir <- sprintf("%s__Param-%s", sub.dir, comb[, "parametric"])
    }

    # If necessary, remove defaults
    if (shortened) {
      components <- str_split(sub.dir, "__")[[1]]
      components.shortened <- NULL
      for (component in components) {
        if (!(str_split(component, "-")[[1]][1] %in% names(defaults))) {
          components.shortened <- c(components.shortened, component)
        }
      }

      sub.dir <- paste(components.shortened, collapse = "__")
    }

    # Create directory if necessary
    check_create.dir(master.dir)
    check_create.dir(sub.dir, master.dir)

    # File name: (gs.method,) (n.if.per.cov,) n.sim, n, n.cov, DGP, B, t.eval
    file.name <- if (comb[, "search.method"] == "GS") {sprintf("method-%s", comb[, "gs.method"])} else NULL
    file.name <- if (comb[, "inst.func.family"] == "box") {
      paste(c(file.name, sprintf("nbpc-%s", comb[, "n.if.per.cov"])), collapse = "__")
    } else if (comb[, "inst.func.family"] == "spline") {
      paste(c(file.name, sprintf("nifpc-%s__degree-%s", comb[, "n.if.per.cov"], comb[, "degree"])), collapse = "__")
    } else if (comb[, "inst.func.family"] %in% c("cd", "cdmc")) {
      paste(c(file.name, sprintf("nifpc-%s__Gc-%s__degree-%s", comb[, "n.if.per.cov"], comb[, "inst.func.family.c"], comb[, "degree"])), collapse = "__")
    } else file.name
    gen.info <- sprintf("geninfo-%s", paste(c(comb[, "n"], comb[, "n.cov"],
                                              comb[, "DGP"], comb[, "B"],
                                              comb[, "t.eval"], seed, sim.iter.nbr),
                                            collapse = "__"))
    file.name <- ifelse(is.null(file.name),
                        gen.info,
                        paste(c(file.name, gen.info), collapse = "__"))
    file.name <- paste0(file.name, ".csv")

    # If necessary, remove defaults
    if (shortened) {
      components <- str_split(file.name, "__")[[1]]
      components.shortened <- NULL
      for (component in components) {
        if (!(str_split(component, "-")[[1]][1] %in% names(defaults))) {
          components.shortened <- c(components.shortened, component)
        }
      }

      file.name <- paste(components.shortened, collapse = "__")
    }

  } else if (is(comb, "list")) {

    # Shorten inst.func.family name, if applicable.
    if (comb[["inst.func.family"]] == "cd.manycov") {
      comb[["inst.func.family"]] <- "cdmc"
    }

    # Directory: search method, inst.func.family, K.bar
    sub.dir <- sprintf("Search-%s__IF-%s__Kbar-%s__Alpha-%s", comb[["search.method"]],
                       comb[["inst.func.family"]], comb[["K.bar"]], comb[["alpha"]])
    if ("link.function" %in% names(comb)) {
      sub.dir <- sprintf("LF-%s__%s", comb[["link.function"]], sub.dir)
    }
    if ("parametric" %in% names(comb)) { # Option used in misspecification simulations
      sub.dir <- sprintf("%s__Param-%s", sub.dir, comb[["parametric"]])
    }

    # If necessary, remove defaults
    if (shortened) {
      components <- str_split(sub.dir, "__")[[1]]
      components.shortened <- NULL
      for (component in components) {
        if (!(str_split(component, "-")[[1]][1] %in% names(defaults))) {
          components.shortened <- c(components.shortened, component)
        }
      }

      sub.dir <- paste(components.shortened, collapse = "__")
    }

    # Create directory if necessary
    check_create.dir(master.dir)
    check_create.dir(sub.dir, master.dir)

    # File name: (gs.method,) (n.if.per.cov,) n.sim, n, n.cov, DGP, B, t.eval
    file.name <- if (comb[["search.method"]] == "GS") {sprintf("method-%s", comb[["gs.method"]])} else NULL
    file.name <- if (comb[["inst.func.family"]] == "box") {
      paste(c(file.name, sprintf("nbpc-%s", comb[["n.if.per.cov"]])), collapse = "__")
    } else if (comb[["inst.func.family"]] == "spline") {
      paste(c(file.name, sprintf("nifpc-%s__degree-%s", comb[["n.if.per.cov"]], comb[["degree"]])), collapse = "__")
    } else if (comb[["inst.func.family"]] %in% c("cd", "cdmc")) {
      paste(c(file.name, sprintf("nifpc-%s__Gc-%s__degree-%s", comb[["n.if.per.cov"]], comb[["inst.func.family.c"]], comb[["degree"]])), collapse = "__")
    } else file.name
    gen.info <- sprintf("geninfo-%s", paste(c(comb[["n"]], comb[["n.cov"]],
                                              comb[["DGP"]], comb[["B"]],
                                              comb[["t.eval"]], seed, sim.iter.nbr),
                                            collapse = "__"))
    file.name <- ifelse(is.null(file.name),
                        gen.info,
                        paste(c(file.name, gen.info), collapse = "__"))
    file.name <- paste0(file.name, ".csv")

    # If necessary, remove defaults
    if (shortened) {
      components <- str_split(file.name, "__")[[1]]
      components.shortened <- NULL
      for (component in components) {
        if (!(str_split(component, "-")[[1]][1] %in% names(defaults))) {
          components.shortened <- c(components.shortened, component)
        }
      }

      file.name <- paste(components.shortened, collapse = "__")
    }
  }

  # Obtain full name
  path <- paste(master.dir, sub.dir, file.name, sep = "/")

  # Return the result
  path
}

#' @title Get matrix of simulation settings
#'
#' @description This function takes a list of simulation parameters as argument
#' and returns a matrix where each row in the matrix corresponds to a setting
#' that should be simulated. This is especially handy if we want to work with
#' the 'worker' module on the VSC. This function is used in 'run.simulations.R',
#' which is not included here.
#'
#' @param sim.params List of simulation paramerers.
#' @param for.VSC Boolean value indicating whether the matrix returned by
#' this function is meant to be used as input to the worker module in the super-
#' computer. If \code{for.VSC = TRUE}, additional columns will be added to the
#' the matrix, the matrix will be saved and will not be returned. Default is
#' \code{for.VSC = FALSE}.
#' @param for.cCDF Boolean value indicating whether the matrix returned by this
#' function is meant to be used in the conditional CDF simulations. If
#' \code{for.cCDF = TRUE}, the matrix will be adapted so that the same seed is
#' used among desings that differ only in the value for \code{t.eval}. Default
#' is \code{for.cCDF = FALSE}.
#'
#' @returns Simulation setting matrix
#'
#' @noRd
#'
get.simulation.variable.matrix <- function(sim.params, for.VSC = FALSE,
                                           for.cCDF = FALSE) {

  # Extract variables used in simulation for readability
  n.sim.vct <- sim.params[["n.sim"]]
  n.vct <- sim.params[["n"]]
  n.cov.vct <- sim.params[["n.cov"]]
  DGP.vct <- sim.params[["DGP"]]
  search.method.vct <- sim.params[["search.method"]]
  alpha.vct <- sim.params[["alpha"]]
  t.eval.vct <- sim.params[["t.eval"]]
  K.bar.vct <- sim.params[["K.bar"]]
  n.if.per.cov.vct <- sim.params[["n.if.per.cov"]]
  B.vct <- sim.params[["B"]]
  inst.func.family.vct <- sim.params[["inst.func.family"]]
  inst.func.family.c.vct <- sim.params[["inst.func.family.c"]]
  degree.vct <- sim.params[["degree"]]
  gs.method.vct <- sim.params[["gs.method"]]
  iseed <- sim.params[["iseed"]]
  idx.param.of.interest <- sim.params[["idx.param.of.interest"]]

  # Extract variables related to the data generation
  beta.true <- sim.params[["beta.true"]]
  par.space <- sim.params[["par.space"]]

  # Master directory for storing simulation results
  master.dir <- sim.params[["master.dir"]]

  # The grid search method might not have been specified when simulations will
  # only take place using the EAM algorithm.
  if (!("GS" %in% search.method.vct)) {
    gs.method.vct <- 1
  }

  # The degree of the spline instrumental functions might not have been specified
  # when the family of instrumental functions to be used are box functions.
  if ((length(inst.func.family.vct) == 1) & ("box" %in% inst.func.family.vct)) {
    degree.vct <- 0
  }

  # The family of instrumental functions to be used when G.cd is selected might
  # not be specified when G.cd is not selected.
  if (!("cd" %in% inst.func.family.vct)) {
    inst.func.family.c.vct <- "ignore"
  }

  # All simulation parameters must be specified
  if (any(is.null(c(n.sim.vct, n.vct, n.cov.vct, DGP.vct, search.method.vct,
                    alpha.vct, gs.method.vct, K.bar.vct, n.if.per.cov.vct,
                    inst.func.family.vct, degree.vct, inst.func.family.c.vct,
                    B.vct, t.eval.vct, iseed, dir, c)))) {
    stop("Some simulation parameters are not specified")
  }

  # If the matrix of combinations is meant to be supplied to the supercomputer...
  if (for.VSC) {

    # Re-code the string-valued options into integers (turns out to be not
    # necessary but I'll just leave it like this).
    search.method2int <- function(sm.str) {
      if (sm.str == "GS") {return(1)}
      if (sm.str == "EAM") {return(2)}
    }
    search.method.vct <- unname(Vectorize(search.method2int)(search.method.vct))
    inst.func.fam2int <- function(iff.str) {
      if (iff.str == "box") {return(1)}
      if (iff.str == "spline") {return(2)}
      if (iff.str == "cd") {return(3)}
      else {return(0)}
    }
    inst.func.family.vct <- unname(Vectorize(inst.func.fam2int)(inst.func.family.vct))
    inst.func.family.c.vct <- unname(Vectorize(inst.func.fam2int)(inst.func.family.c.vct))

    # For simplicity, only allow for one option for the number of simulations to
    # be run.
    if (length(n.sim.vct) != 1) {
      stop("n.sim.vct must have length 1.")
    }
    n.sim <- n.sim.vct

    # Obtain matrix of all parameter combinations to run
    combinations <- expand.grid(n.vct, n.cov.vct, DGP.vct, search.method.vct,
                                alpha.vct, gs.method.vct, K.bar.vct, n.if.per.cov.vct,
                                inst.func.family.vct, degree.vct, inst.func.family.c.vct,
                                B.vct, t.eval.vct)
    colnames(combinations) <- c("n", "n.cov", "DGP", "search.method", "alpha",
                                "gs.method", "K.bar", "n.if.per.cov",
                                "inst.func.family", "degree", "inst.func.family.c",
                                "B", "t.eval")

    # Initialize vector of seeds to be used
    seeds.to.use <- iseed + 1:n.sim

    # Each simulation setting should be run n.sim amount of times, all with a
    # different initial seed.
    comb.extended <- matrix(nrow = 0, ncol = ncol(combinations) + 1)
    colnames(comb.extended) <- c("seed", colnames(combinations))
    comb.extended <- as.data.frame((comb.extended))
    for (i in 1:nrow(combinations)) {
      for (seed in seeds.to.use) {
        row <- c(seed = seed, combinations[i, ])
        comb.extended <- rbind(comb.extended, row)
      }
      seeds.to.use <- seeds.to.use + n.sim
    }

    # When the search method is not 'grid search', we should not simulate over
    # different values for 'grid search method'.
    comb.extended[which(comb.extended$search.method != "GS"), "gs.method"] <- 1

    # When the family of instrumental functions is the box family, we should not
    # simulate over different degrees.
    comb.extended[which(comb.extended$inst.func.family == 1), "degree"] <- 0

    # When the family of instrumental functions is not the continuous/discrete
    # family, we should not iterate over different families for the continuous
    # part.
    comb.extended[which(comb.extended$inst.func.family != 3), "inst.func.family.c"] <- 0

    # When the family of instrumental functions is the continuous/discrete
    # family and the family for the continuous part are box functions, we should
    # not iterate over different values of degree.
    comb.extended[which((comb.extended$inst.func.family == 3) & (comb.extended$inst.func.family.c == 1)), "degree"] <- 0

    # When the degree of the instrumental functions is not smaller than the
    # number of instrumental functions to be used per covariate, the method will
    # throw an error. Remove such cases.
    comb.extended <- comb.extended[!(comb.extended$n.if.per.cov <= comb.extended$degree), ]

    # Given the previous modifications, obtain all unique simulation settings.
    comb.extended <- unique(comb.extended)

    # Rename columns
    colnames(comb.extended) <- c("seed", "n", "n_cov", "DGP", "search_method", "alpha",
                                 "gs_method", "K_bar", "n_if_per_cov",
                                 "inst_func_family", "degree", "inst_func_family_c",
                                 "B", "t_eval")

    # If necessary, adapt the matrix for use in cCDF simulations
    if (for.cCDF) {

      # Get all simulation designs: remove seed and t_eval information
      designs <- comb.extended[, !(colnames(comb.extended) %in% c("seed", "t_eval"))]

      # Get all unique simulation designs
      unique.designs <- unique(designs)

      # For each unique simulation design, set the same seed to be used over all
      # values of t.eval.
      for (design.idx in 1:nrow(unique.designs)) {

        # Select unique desing of this iteration
        design <- unique.designs[design.idx,]

        # Get row indices of all simulations corresponding to this design
        equal.design.idxs <- which(apply(designs, 1, function(row) {all(row == design)}))

        # For each, make the seed to be used equal
        seed.to.use <- comb.extended[min(equal.design.idxs), "seed"]
        comb.extended[equal.design.idxs, "seed"] <- seed.to.use
      }
    }

    # Save the results
    write.csv(comb.extended, "combinations.csv", row.names = FALSE)
  }

  # If the matrix of combinations is not meant to be supplied to the super-
  # computer...
  if (!for.VSC) {

    # Obtain matrix of all parameter combinations to run
    combinations <- expand.grid(n.sim.vct, n.vct, n.cov.vct, DGP.vct, search.method.vct,
                                alpha.vct, gs.method.vct, K.bar.vct, n.if.per.cov.vct,
                                inst.func.family.vct, degree.vct, inst.func.family.c.vct,
                                B.vct, t.eval.vct, stringsAsFactors = FALSE)
    colnames(combinations) <- c("n.sim", "n", "n.cov", "DGP", "search.method",
                                "alpha", "gs.method", "K.bar", "n.if.per.cov",
                                "inst.func.family", "degree", "inst.func.family.c",
                                "B", "t.eval")

    # When the search method is not 'grid search', we should not simulate over
    # different values for 'grid search method'.
    combinations[which(combinations$search.method != "GS"), "gs.method"] <- 1

    # When the family of instrumental functions is the box family, we should not
    # simulate over different degrees.
    combinations[which(combinations$inst.func.family == "box"), "degree"] <- 0

    # When the family of instrumental functions is not the continuous/discrete
    # family, we should not iterate over different families for the continuous
    # part.
    combinations[which(combinations$inst.func.family != "cd"), "inst.func.family.c"] <- "NA"

    # When the family of instrumental functions is the continuous/discrete
    # family and the family for the continuous part are box functions, we should
    # not iterate over different values of degree.
    combinations[which((combinations$inst.func.family == "cd") & (combinations$inst.func.family.c == "box")), "degree"] <- 0

    # Given the previous modifications, obtain all unique simulation settings.
    combinations <- unique(combinations)

    # Add the initial seeds to be used
    iseed.col <- rep(iseed, nrow(combinations))
    n.sim <- 0
    for (i in 1:nrow(combinations)) {
      iseed <- iseed + n.sim
      iseed.col[i] <- iseed
      n.sim <- combinations[i, "n.sim"]
    }
    iseed <- iseed.col
    combinations <- cbind(combinations, iseed)

    # If necessary, adapt the matrix for use in cCDF simulations
    if (for.cCDF) {

      # Get all simulation designs: remove seed and t_eval information
      designs <- combinations[, !(colnames(combinations) %in% c("iseed", "t.eval"))]

      # Get all unique simulation designs
      unique.designs <- unique(designs)

      # For each unique simulation design, set the same seed to be used over all
      # values of t.eval.
      for (design.idx in 1:nrow(unique.designs)) {

        # Select unique desing of this iteration
        design <- unique.designs[design.idx,]

        # Get row indices of all simulations corresponding to this design
        equal.design.idxs <- which(apply(designs, 1, function(row) {all(row == design)}))

        # For each, make the seed to be used equal
        seed.to.use <- combinations[min(equal.design.idxs), "iseed"]
        combinations[equal.design.idxs, "iseed"] <- seed.to.use
      }
    }

    # Return the results
    return(combinations)
  }
}

#' @title Run a simulation for a specific set of variables
#'
#' @description This function is used in 'run.simulations.R', which is not
#' included here.
#'
#' @param comb Data frame of one single row containing the parameter settings to
#' be used in the simulation.
#' @param beta.true True values of the parameters in the model; used to generate
#' the data.
#' @param idx.param.of.interest Index of the element of the parameter vector one
#' wants to create the identified interval for. I.e. if c is the projection
#' vector, then c[idx.param.of.interest] = 1, and c = 0 elsewhere.
#' @param par.space Bounds on the parameter space.
#' @param starting.seed Initial random seed to use
#' @param master.dir Name of the directory in which the simulation results
#' should be stored.
#' @param verbose Verbosity parameter.
#'
#' @importFrom utils write.csv
#'
#' @note [ToDo] Also allow for simulation using other link functions
#'
#' @noRd
#'
simulate1D <- function(comb, beta.true, idx.param.of.interest, par.space,
                       starting.seed, master.dir, verbose) {

  #### Extract hyperparameters ####

  # Extract the parameters (currently, cov.ranges will never be in 'comb')
  n.sim <- comb[, "n.sim"]
  n <- comb[, "n"]
  n.cov <- comb[, "n.cov"]
  DGP <- comb[, "DGP"]
  search.method <- comb[, "search.method"]
  alpha <- comb[, "alpha"]
  t.eval <- comb[, "t.eval"]
  K.bar <- comb[, "K.bar"]
  n.if.per.cov <- comb[, "n.if.per.cov"]
  B <- comb[, "B"]
  inst.func.family <- comb[, "inst.func.family"]
  degree <- comb[, "degree"]
  inst.func.family.c <- comb[, "inst.func.family.c"]
  gs.method <- comb[, "gs.method"]
  cov.ranges <- if ("cov.ranges" %in% colnames(comb)) {comb[, "cov.ranges"]} else NULL

  # Precondition checks. If at some point this would be adapted, also adapt the
  # code in 'conditionalCDFFunctions.R -> get.Peterson.bounds.R'
  beta0.checks <- simplify2array(lapply((-10):10, function(x) {beta.true(x)[1]}))
  if (any(beta0.checks != (-10):10)) {
    stop("Currently only able to handle beta0 = t. (*)")
    # Technically, beta0 = t + b for any b would also work.
  }

  # Projection vector that projects the parameter vector onto the element of
  # interest
  c <- rep(0, n.cov + 1)
  c[idx.param.of.interest] <- 1

  # Family of instrumental functions to be used when G.cd is selected
  G.c <- NULL
  if (inst.func.family == "cd") {
    if (inst.func.family.c == "box") {
      G.c <- G.box
    } else if (inst.func.family.c == "spline") {
      G.c <- G.spline
    } else {
      stop("Specified family of instrumental functions not implemented.")
    }
  }

  # If applicable, select the type of grid search to be carried out
  if (gs.method == 1) {
    next.gs.point = gs.binary  # Binary search
  } else if (gs.method == 2) {
    next.gs.point = gs.regular # Regular grid search
  }

  # List of hyperparameters to be supplied to various functions later on
  options <- list(n.if.per.cov = n.if.per.cov,
                  K.bar = K.bar,
                  B = B,
                  gs.method = gs.method,
                  next.gs.point = next.gs.point,
                  alpha = alpha,
                  link.function = "AFT_ll",
                  DGP = DGP,
                  inst.func.family = inst.func.family,
                  degree = degree,
                  G.c = G.c,
                  cov.ranges = cov.ranges)

  # Flag for recording the time that each step of the algorithm takes
  time.run.duration <- TRUE

  #### Perform the simulations ####

  # Initialize object that will store the results
  ncol.sim.results <- 7 + 3 * as.numeric(search.method == "EAM") +
    2 * as.numeric(DGP %% 20 %in% c(10, 11, 12))
  sim.results <- matrix(nrow = 0, ncol = ncol.sim.results)
  cols <- c("ident.set.l", "ident.set.u", "conv.l", "conv.u", "total.run.time")
  if (search.method == "EAM") {
    cols <- c(cols, "E-step", "A-step", "M-step")
  }
  cols <- c(cols, "seed", "per.cens")
  if ((DGP %% 20) %in% c(10, 11, 12)) {
    cols <- c(cols, "per.cens.reg1", "per.cens.reg2")
  }
  colnames(sim.results) <- cols

  # Run the simulation
  for (sim.iter.nbr in 1:n.sim) {
    if (verbose >= 1) {
      message(sprintf("Starting simulation %d / %d", sim.iter.nbr, n.sim))
    }

    # Set random seed of this iteration
    seed.to.use <- starting.seed + sim.iter.nbr
    set.seed(seed.to.use)

    # Generate data
    data <- generateData(beta.true, n, n.cov, options, plot_data = (verbose >= 3))

    # Find the identified set based on Bei (2024)
    fis.out <- find.identified.set(c, t.eval, par.space, data, search.method,
                                   options, verbose, time.run.duration)
    ident.set <- fis.out$ident.set

    # Analyse the running times
    chrono1 <- fis.out$chronometer1
    chrono2 <- fis.out$chronometer2
    total.run.time <- chrono1$get.total.time(force = TRUE) + chrono2$get.total.time(force = TRUE)
    if (search.method == "EAM") {
      total.leg.times <- chrono1$accumulate.legs(force = TRUE)[1:3] +
        chrono2$accumulate.legs(force = TRUE)[1:3]
    } else {
      total.leg.times <- NULL
    }

    # Convergence of the algorithm
    conv.l <- fis.out$converge2
    conv.u <- fis.out$converge1

    # Censoring percentage in the data
    per.cens <- 1 - sum(data$Delta)/n
    per.cens.region1 <- NULL
    per.cens.region2 <- NULL
    if ((DGP %% 20) %in% c(10, 11, 12)) {
      per.cens.region1 <- 1 - sum(data$Delta[inRegion1(data)])/length(inRegion1(data))
      per.cens.region2 <- 1 - sum(data$Delta[inRegion2(data)])/length(inRegion2(data))
    }

    # Store the results
    sim.result <- c(ident.set, conv.l, conv.u, total.run.time, total.leg.times,
                    seed.to.use, per.cens, per.cens.region1, per.cens.region2)
    sim.results <- rbind(sim.results, sim.result)

    # Store (also intermediate) results
    path <- get.file.name(comb, sim.iter.nbr, seed.to.use, master.dir)
    write.csv(sim.results, path, row.names = FALSE)
  }

  # Return the number of simulations ran (used in updating the initial seed for
  # the next simulation)
  n.sim
}

#' @title Estimate an identified model assuming independence.
#'
#' @description This function estimates an identified model assuming that the
#' event and censoring time independent. Specifically, it estimates a regular
#' AFT model or Cox PH model. This is used as a reference in the simulations.
#'
#' @param c Projection vector
#' @param t.eval Time point of interest.
#' @param data Data frame.
#' @param link.function Name of the link function to use ("AFT_ll" or "Cox_wb").
#'
#' @returns The coefficients of the identified model.
#'
#' @importFrom survival coxph survreg
#' @importFrom stats as.formula
#'
#' @noRd
#'
get.coef.identified.model <- function(c, t.eval, data, link.function) {

  # Get number of covariates
  n.cov <- sum(grepl("X[1-9][[:digit:]]*$", colnames(data)))

  # Estimate the model
  if (link.function == "AFT_ll") {
    fml <- as.formula(paste0("Surv(Y, Delta) ~ ", paste(paste0("X", 1:n.cov), collapse = " + ")))
    fit <- survreg(fml, dist = "loglogistic", data = data)
    coef <- (-fit$coefficients) %*% c

  } else if (link.function == "Cox_wb") {
    fml <- as.formula(paste0("Surv(Y, Delta) ~ ", paste(paste0("X", 1:n.cov), collapse = " + ")))
    fit <- coxph(fml, data = data)
    coef <- fit$coefficients %*% c[-1]
  }

  # Return the coefficient of interest
  coef
}



