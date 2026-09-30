
# This file contains an R translation of the .ado code file that can be found
# on 'https://sites.google.com/view/brighamfrandsen/software#h.lhxa3mg5af8s'.

#' @title Compute the reduced-sample estimator.
#'
#' @description
#' This function computes F^{RS}(y|x) as shown in Equation (3) of Frandsen (2019).
#'
#' @param y Value at which to compute the estimator.
#' @param Y Vector of observed times.
#' @param C Vector of censoring times.
#' @param cell.vec full vector of cell indices, for each observation.
#' @param cell Cell for which the estimator should be computed (corresponding
#' to the value of x on which is conditioned in the expression F^{RS}(y|x)).
#'
#' @noRd
#'
F.RS <- function(y, Y, C, cell.vec, cell) {

  # Compute numerator and denominator
  n <- sum(as.numeric(Y <= y) * as.numeric((C > y) & (cell.vec == cell)))
  d <- sum(as.numeric((C > y) & (cell.vec == cell)))

  # Return the result
  n/(d + 1e-10)
}

#' @title Compute the Kaplan-Meier weights.
#'
#' @description
#' This function computes the formula as defined in Step 1 of the algorithm
#' outline provided by Frandsen (2019), Section 3.2 - Test Procedure.
#'
#' @noRd
#'
W.KM <- function(i, Y, C, nx) {
  event <- as.numeric(Y < C)
  denom <- nx - i + 1
  j.vec <- 1:(i - 1)
  if (i >= 2) {
    fac <- prod( ( (nx - j.vec) / (nx - j.vec + 1) )^event[1:(i - 1)])
  } else {
    fac <- 1
  }
  event[i] / denom * fac
}

#' @title Compute the test statistic of Frandsen (2019).
#'
#' @description
#' This function computes the test statistic displayed in Equation (4) of
#' Frandsen (2019).
#'
#' @param data Data frame.
#' @param Y.name Name of the variable in \code{data} that represents the
#' observed time min(T, C).
#' @param C.name Name of the variable in \code{data} that represents the
#' censoring times.
#' @param X.name Name of the discrete variable in \code{data} that represents
#' the single discrete covariate.
#' @param TRIM_POINT Quantile of C at which to trim the integral when computing
#' the test statistic. Frandsen recommends the value 0.75, which is therefore
#' selected as the default value.
#' @param warn Boolean value indicating whether warnings should be communicated
#' to the user. Default is \code{warn = TRUE}.
#'
#' @importFrom stats quantile
#'
#' @noRd
#'
calcstat <- function(data, Y.name, C.name, X.name = NULL,
                     TRIM_POINT = 0.75, warn = TRUE) {

  # Require the survival package
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required to run the test of Frandsen (2019).")
  }

  # Divide the data into cells based on the covariate value. If no covariate
  # is provided, all data will be placed into cell 1.
  if (is.null(X.name)) {
    data$cell <- 1
  } else {
    data$cell <- interaction(data[, X.name], drop = TRUE)
  }
  cells <- unique(data$cell)

  # Inform users if the number of cells is large (indicative of the presence of
  # a continuous covariate)
  if (length(cells) > max(10, sqrt(nrow(data))) & warn) {
    warning("Large number of unique covariate values detected. ",
            "Note that only discrete covariates may be supplied!")
  }

  # Define some useful variables (sample size; trim point -> See Frandsen, 2019)
  n <- nrow(data)
  tau <- TRIM_POINT

  # Compute the test statistic Delta (Equation (4) in Frandsen).
  delta <- 0
  for(cell in cells) {

    # Subset data to cell of this iteration.
    data.cell <- data[data$cell == cell,]

    # Compute sample size of data subset.
    nx <- nrow(data.cell)

    # Obtain the cut-off time
    maxt <- quantile(data.cell[[C.name]], tau)

    # Reduced-sample estimator of distribution function
    Fc <- unlist(lapply(data.cell[[Y.name]], F.RS, Y = data[[Y.name]],
                        C = data[[C.name]], cell.vec = data$cell, cell = cell))

    # Kaplan-Meier estimator of distribution function
    event <- as.numeric(data.cell[[Y.name]] != data.cell[[C.name]])
    fit <- survival::survfit(survival::Surv(data.cell[[Y.name]], event) ~ 1)
    KM <- 1 - summary(fit, times = data.cell[[Y.name]], extend = TRUE)$surv

    # Sort both by Y
    Y.order <- order(data.cell[[Y.name]])
    Y <- data[[Y.name]][Y.order]
    C <- data[[C.name]][Y.order]
    Fc <- Fc[Y.order]
    KM <- KM[Y.order]

    # Compute the difference
    diffs <- Fc - KM

    # Compute the Kaplan-Meier weights
    KM.weights <- unlist(lapply(1:nx, W.KM, Y = Y, C = C, nx = nx))

    # Compute the test statistic for this cell
    T.stat.cell <- sum(diffs^2 * KM.weights)

    # Add the test statistic computed on the cell to the full test statistsic,
    # weighted by the probability of X beloning to that cell.
    delta <- delta + (nx/n) * T.stat.cell
  }

  # Return the result.
  delta
}

#' @title Perform test of Frandsen (2019).
#'
#' @description
#' This function computes the test of Frandsen (2019). It is based on
#' subsampling and has very much unoptimized implementation, so can be rahter
#' slow.
#'
#' @param data Data frame.
#' @param Y.name Name of the variable in \code{data} that represents the
#' observed time min(T, C).
#' @param C.name Name of the variable in \code{data} that represents the
#' censoring times.
#' @param X.name Name of the discrete variable in \code{data} that represents
#' the single discrete covariate allowed in the testing procedure.
#' @param ssreps Number of subsampling repetitions to perform.
#' @param TRIM_POINT Quantile of C at which to trim the integral when computing
#' the test statistic. Frandsen recommends the value 0.75, which is therefore
#' selected as the default value.
#'
#' @references Brigham R. Frandsen (2019) Testing Censoring Point Independence, Journal
#' of Business & Economic Statistics, 37:3, 496-505, DOI: 10.1080/07350015.2017.1383261
#'
#' @examples
#' \donttest{
#' # Generate survival data with censoring time always observed
#' n <- 5000
#' X <- sample(0:1, 10, replace = TRUE)
#' U <- copula::rCopula(n, copula::frankCopula(param = 6))
#' T <- X + qexp(U[, 1], rate = 1)
#' C <- X + qexp(U[, 2], rate = 1.5)
#' Y <- pmin(T, C)
#' Delta <- as.numeric(Y == T)
#' data <- as.data.frame(cbind(Y, Delta, C, X))
#' colnames(data) <- c("Y", "Delta", "C", "X")
#'
#' # Run test of Frandsen
#' Y.name <- "Y"
#' C.name <- "C"
#' X.name <- "X"
#' ssreps <- 500
#' TRIM_POINT <- 0.75
#' FrandsenTest(data, Y.name, C.name, X.name, ssreps, TRIM_POINT)
#' }
#'
#' @export
#'
FrandsenTest <- function(data, Y.name, C.name, X.name = NULL,
                         ssreps = 500, TRIM_POINT = 0.75) {

  # Require the survival package
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required to run the test of Frandsen (2019).")
  }

  # Define variable names
  vars <- c(Y.name, C.name, X.name)

  # Check for missing values. If there are any, throw an error
  if (sum(complete.cases(data[, vars])) != nrow(data)) {
    stop("Missing values in relevant columns found.")
  }

  # Subset data to necessary variables.
  data <- data[complete.cases(data[, vars]), vars]
  fullsamplesize <- nrow(data)

  # full sample statistic
  delta <- calcstat(data, Y.name, C.name, X.name, TRIM_POINT = TRIM_POINT)

  # subsample size
  minfeas <- 20
  B <- ceiling(minfeas + fullsamplesize^0.75)
  numsubs <- min(choose(fullsamplesize,B), ssreps)

  # For each subsample iteration, compute the test statistic
  ssdelta <- numeric(numsubs)
  for(b in 1:numsubs) {

    # Get subsample
    pick <- sample(1:fullsamplesize, B)
    subsample <- data[pick,]

    # Compute test statistic on subsample
    stat <- calcstat(subsample, Y.name, C.name, X.name, TRIM_POINT = TRIM_POINT,
                     warn = FALSE)

    # Store result
    ssdelta[b] <- stat
  }

  # Compute test statistic
  ssdelta <- B * (ssdelta - delta)
  delta_scaled <- fullsamplesize * delta
  pvalue <- mean(round(ssdelta,7) >= round(delta_scaled,7), na.rm = TRUE)

  # Return the result
  list(
    test_statisic = delta_scaled,
    p_value = pvalue,
    sample_size = fullsamplesize
  )
}
