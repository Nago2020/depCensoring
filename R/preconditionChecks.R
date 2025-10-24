
#' @title Standardize data format
#'
#' @description Checks the required preconditions of the data and possibly
#' restructures the data.
#'
#' @param data A data frame that should contain columns named \code{Y} and
#' \code{delta} (unless \code{comp.risks = TRUE}, see later).
#' @param admin Boolean value indicating whether the provided data frame contains
#' administrative (i.e. independent) censoring on top of the dependent censoring
#' (in the column named \code{delta}). The default is \code{admin = FALSE}.
#' @param conf Boolean value indicating whether the provided data frame contains
#' a confounded variable and a corresponding instrument. If \code{cond = TRUE},
#' the provided data frame should contain columns named \code{Z} and \code{W},
#' corresponding to the confounded variable and instrument, respectively.
#' Moreover, \code{Zbin} and \code{Wbin} should be specified. The default value
#' is \code{conf = FALSE}.
#' @param comp.risks Boolean value indicating whether the provided data frame
#' contains competing risks. If \code{comp.risks = TRUE}, the given data frame
#' should contain the columns \code{delta1}, \code{delta2}, etc., corresponding
#' to the indicators I(Y = T1), I(Y = T2), etc. respectively. The default is
#' \code{comp.risks = FALSE}.
#' @param Zbin Boolean or integer value (0, 1) indicating whether the confounded
#' variable is binary. \code{Zbin = TRUE} or \code{Zbin = 1} means that Z is
#' binary. \code{Zbin = FALSE} or \code{Zbin = 0} means that Z is continuous.
#' @param Wbin Boolean or integer value (0, 1) indicating whether the instrument
#' is binary. \code{Wbin = TRUE} or \code{Wbin = 1} means that W is binary.
#' \code{Wbin = FALSE} or \code{Wbin = 0} means that W is continuous.
#'
#' @return Returns the uniformized data set.
#'
#' @noRd
#'


uniformize.data <- function(data, admin = FALSE, conf = FALSE, comp.risks = FALSE,
                            Zbin = NULL, Wbin = NULL) {

  # ToDo?: Option to turn off precondition checking .

  #### Step 1: Precondition checks ####

  # Put all column names to lower case
  colnames(data) <- tolower(colnames(data))

  # Given data object should be a data frame
  if (!is.data.frame(data)) {
    stop("The data should be provided as a data frame")
  }

  ## Data cannot already contain an intercept

  if (any(apply(data == 1, MARGIN = 2, all))) {
    stop("The given data frame cannot already contain an intercept")
  }

  ## Necessary column names

  # The time variable
  if (!("y" %in% colnames(data))) {
    stop("The column corresponding to the time variable should be named 'Y'.")
  }

  # Censoring indicators
  if (comp.risks) {
    if ("delta" %in% colnames(data)) {
      stop(paste0("Column named 'delta' detected in data frame. Since comp.risks",
                  " = TRUE, columns containing censoring indicators should be ",
                  "named 'deltak', where k = 1, 2, ..."))
    }
    if ("xi" %in% colnames(data)) {
      stop(paste0("Column named 'xi' detected in data frame. Since comp.risks",
                  " = TRUE, naming conventions for censoring indicator columns ",
                  "are different."))
    }
    if (!("delta1" %in% colnames(data))) {
      stop(paste0("The column corresponding to the indicator I(Y = T1) should ",
                  "be named 'delta1' and analogous naming should be used for ",
                  "the column corresponding to I(Y = Ti). No column 'delta1' ",
                  "detected."))
    }
    if (admin & !("da" %in% colnames(data))) {
      stop(paste0("The column corresponding to the indicator I(Y = A) should be ",
                  "named 'da'. If there is no administrative censoring in the ",
                  "data (and hence this column doesn't exist), please set ",
                  "admin = FALSE."))
    }
  } else {
    if (!("delta" %in% colnames(data))) {
      stop("The column corresponding to the indicator I(Y = T) should be named 'delta'.")
    }
    if (admin & !("xi" %in% colnames(data))) {
      stop("The column corresponding to indicator I(Y = C) should be named 'xi'.")
    }
    if (any(grepl('delta[[:digit:]]+', colnames(data)))) {
      stop(paste0("Column(s) named 'deltak' detected for some integer(s) k. ",
                  "Since comp.risks = FALSE, naming conventions for censoring ",
                  "indicators are different."))
    }
    if ("da" %in% colnames(data)) {
      stop(paste0("Column named 'da' detected. ",
                  "Since comp.risks = FALSE, naming conventions for censoring ",
                  "indicators are different."))
    }
  }

  if (conf & !all(c("z", "w") %in% colnames(data))) {
    stop(paste0("The column corresponding to the confounded variable should be",
                " named 'Z' and the column corresponding to the instrument",
                " should be named 'W'."))
  }

  # Indices of columns corresponding to covariates
  if (comp.risks) {
    nbr.comp.risks <- length(which(substr(colnames(data), 1, 5) == "delta"))
    cov.col.idxs <- setdiff(1:ncol(data), which(colnames(data) %in% c("y", "delta", "xi", paste0("delta", 1:nbr.comp.risks), "da", "z", "w")))
  } else {
    cov.col.idxs <- setdiff(1:ncol(data), which(colnames(data) %in% c("y", "delta", "xi", "z", "w")))
  }


  ## Consistent arguments
  if (!is.numeric(data[, "y"])) {
    stop("Column corresponding to 'Y' should be numeric")
  }
  for (cens.col.idx in which(substr(colnames(data), 1, 5) == "delta")) {
    if (!all(data[, cens.col.idx] %in% c(0, 1))) {
      stop(paste0("Censoring indicator in column ", cens.col.idx, " can only ",
                  "take values 0 or 1."))
    }
  }
  if (!all(sapply(data, class)[cov.col.idxs] %in% c("numeric", "integer"))) {
    stop("The covariates should be numeric")
  }

  if (conf) {
    if (any(is.null(Zbin), is.null(Wbin))) {
      stop("When conf == TRUE, the arguments Zbin and Wbin should be provided.")
    }
    if (!((is.numeric(Zbin) | is.logical(Zbin)) & (is.numeric(Wbin) | is.logical(Wbin)))) {
      stop("Zbin and Wbin should be logical or binary")
    }
    if (is.numeric(Zbin) & !(Zbin %in% c(0, 1))) {
      stop("Invalid value for Zbin. Should be either 0 or 1.")
    }
    if (is.numeric(Wbin) & !(Wbin %in% c(0, 1))) {
      stop("Invalid value for Wbin Should be either 0 or 1.")
    }
    if (!all(is.numeric(data[, "z"]), is.numeric(data[, "w"]))) {
      stop("Confounded variable and/or instrument should be numeric")
    }
    if (Zbin & !all(data[, "z"] %in% c(0, 1))) {
      stop("Confounded variable is not a binary as indicated")
    }
    if (Wbin & !all(data[, "w"] %in% c(0, 1))) {
      stop("Instrumental variable is not a binary as indicated")
    }
  } else {
    if (!is.null(c(Zbin, Wbin))) {
      warning("Argument(s) Zbin and/or Wbin ignored as conf == FALSE")
      Zbin <- NULL
      Wbin <- NULL
    }
  }

  if (comp.risks) {
    cens.col.idxs <- which(colnames(data) %in% paste0("delta", 1:nbr.comp.risks))
    if (admin) {
      cens.col.idxs <- c(cens.col.idxs, which(colnames(data) == "da"))
      if (!all(data[, "da"] %in% c(0, 1))) {
        stop("Administrative censoring indicator can only take values 0 or 1.")
      }
    }
    if (!all(apply(X = data[, cens.col.idxs], MARGIN = 1, FUN = sum) == 1)) {
      stop("Precisely one censoring indicator should be 1 for each observation.")
    }
  } else {
    if (admin) {
      if (!all(data[, "xi"] %in% c(0, 1))) {
        stop("Administrative censoring indicator can only take values 0 or 1.")
      }
      if (!all((data[, "delta"] + data[, "xi"]) %in% c(0, 1))) {
        stop("Observations cannot be both administratively and dependently censored at once")
      }
    }
  }

  # Warn for possible user mistakes
  if ((!conf) & any(c("z", "w") %in% colnames(data))) {
    warning(paste0("Possible column(s) corresponding to a confounded variable",
                   " and/or instrument detected in data but conf == FALSE. ",
                   "Treated as exogenous covariates."))
  }
  if ((!conf) & any(c(!is.null(Zbin), !is.null(Wbin)))) {
    warning(paste0("Arguments Zbin and/or Wbin ignored since conf == FALSE"))
  }
  if (conf & !is.null(Zbin)) {
    if (!Zbin) {
      if (all(data[, "z"] %in% c(0, 1))) {
        warning("Z specified as numeric but appears to be binary")
      }
    }
  }
  if (conf & !is.null(Wbin)) {
    if (!Wbin) {
      if (all(data[, "w"] %in% c(0, 1))) {
        warning("W specified as numeric but appears to be binary")
      }
    }
  }
  if (comp.risks) {
    if ((!admin) & ("da" %in% colnames(data))) {
      warning(paste0("Possible administrative censoring column detected in data",
                     " but admin == FALSE so not taken into account"))
    }
  } else {
    if ((!admin) & ("xi" %in% colnames(data))) {
      warning(paste0("Possible administrative censoring column detected in data",
                     " but admin == FALSE so not taken into account"))
    }
  }

  # Missing values
  if (any(is.na(data))) {
    stop("Missing values detected")
  }

  #### Step 2: Create uniformized data frame ####

  # Add intercept column to the data

  data[, "intercept"] <- rep(1, nrow(data))

  if (comp.risks) {

    # Column index of the observed time
    unif.col.idxs <- which(colnames(data) == "y")

    # Column indices of the censoring indicators.
    ind.col.idx <- which(substr(colnames(data), 1, 5) == "delta")

    # Append column indices of censoring indicators to main vector in the order
    # delta1, delta2, delta3, etc.
    unif.col.idxs <- c(unif.col.idxs, ind.col.idx[order(colnames(data)[ind.col.idx])])

    # If administrative censoring: append the column index of indicator I(Y = A)
    # to the main vector.
    if (admin) {
      unif.col.idxs <- c(unif.col.idxs, which(colnames(data) == "da"))
    }

    # Append intercept column to the data
    unif.col.idxs <- c(unif.col.idxs, ncol(data))

    # Append exogenous covariate column indices to main vector.
    unif.col.idxs <- c(unif.col.idxs, cov.col.idxs)

    # If confounding: append the endogenous covariate Z and instrument column W
    # indices to main vector.
    if (conf) {
      unif.col.idxs <- c(unif.col.idxs,  which(colnames(data) == "z"),  which(colnames(data) == "w"))
    }

    # Reorder the data set in specified order.
    unif.data <- data[, unif.col.idxs]
    colnames(unif.data) <- c("y",
                             sort(colnames(data)[ind.col.idx]),
                             if (admin) "da" else NULL,
                             paste0("x", 0:length(cov.col.idxs)),
                             if (conf) "z" else NULL,
                             if (conf) "w" else NULL)

  } else {

    # Column indices of observed time Y and indicator I(Y = T)
    unif.col.idxs <- c(which(colnames(data) == "y"), which(colnames(data) == "delta"))

    # If administrative censoring: append the column index of indicator I(Y = A)
    # to the main vector.
    if (admin) {
      unif.col.idxs <- c(unif.col.idxs, which(colnames(data) == "xi"))
    }

    # Append intercept column to the data
    unif.col.idxs <- c(unif.col.idxs, ncol(data))

    # Append exogenous covariate column indices to main vector.
    unif.col.idxs <- c(unif.col.idxs, cov.col.idxs)

    # If confounding: append the endogenous covariate Z and instrument column W
    # indices to main vector.
    if (conf) {
      unif.col.idxs <- c(unif.col.idxs,  which(colnames(data) == "z"),  which(colnames(data) == "w"))
    }

    # Reorder the data set in specified order.
    unif.data <- data[, unif.col.idxs]
    colnames(unif.data) <- c("y",
                             "delta",
                             if (admin) "xi" else NULL,
                             paste0("x", 0:length(cov.col.idxs)),
                             if (conf) "z" else NULL,
                             if (conf) "z" else NULL)
  }

  return(unif.data)

}
