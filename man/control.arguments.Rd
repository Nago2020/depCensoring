% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SemiParametricCopulawithoutCovariates.R
\name{control.arguments}
\alias{control.arguments}
\title{Prepare initial values within the control arguments}
\usage{
control.arguments(maxit = 300, eps = 1e-06, trace = TRUE, ktau.inits = NULL)
}
\arguments{
\item{maxit}{a positive integer that denotes the maximum iteration number in optimization.}

\item{eps}{a positive small numeric value that denotes the tolerance for convergence.}

\item{trace}{a logical value that judges whereh the tracing information on the progress of the model fitting should be produced. The default value if \code{TRUE}.}

\item{ktau.inits}{a numeric vector that contains initial values of the Kendall's tau.}
}
\description{
Prepare initial values within the control arguments
}
