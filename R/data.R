# Documentation for data sets.

#' @title Liver cirrhosis data set.
#'
#' @description End stage liver disease data set, as for example analyzed by
#' D'Haen et al. (2025).
#'
#' @format A data frame with 281 rows and 7 variables:
#' \describe{
#'  \item{patient}{patient ID.}
#'  \item{time}{Survival time.}
#'  \item{status}{Censoring indicator.}
#'  \item{age}{Age of patient.}
#'  \item{gender}{Gender of patient (1 - male, 0 - female).}
#'  \item{bmi}{Body mass index.}
#'  \item{ukeld}{UK End-stage Liver Disease score. Higher scores indicate more
#'  severe disease state.}
#' }
"liver"
