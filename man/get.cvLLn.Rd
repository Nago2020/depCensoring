% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoundingCovariateEffects.R
\name{get.cvLLn}
\alias{get.cvLLn}
\title{Compute the critical value of the test statistic.}
\usage{
get.cvLLn(
  BetaI.r,
  data,
  t,
  hp,
  c,
  r,
  par.space,
  inst.func.evals = NULL,
  alpha = 0.95
)
}
\arguments{
\item{BetaI.r}{Matrix containing in its columns the minimizers of the
S-function leading to the test statistic.}

\item{data}{Data frame.}

\item{t}{Time point of interest. Also allowed to
be a vector of time points (used in estimating the model under assumed time-
independent coefficients).}

\item{hp}{List of hyperparameters.}

\item{c}{Projection vector.}

\item{r}{Result of projection of parameter vector onto \code{c}.}

\item{par.space}{Bounds on the parameter space.}

\item{inst.func.evals}{Matrix of precomputed instrumental function
evaluations for each observation in the data set. If \code{NULL}, the
evaluations will be computed during execution of this function. Default is
\code{inst.func.evals = NULL}.}

\item{alpha}{Confidence level.}
}
\value{
The critical value for the test statistic.
}
\description{
This function computes the critical value following the
algorithm of Section 4.3 in Bei (2024).
}
\references{
Bei, X. (2024). Local linearieation based subvector inference in
moment inequality models. Journal of Econometrics. 238:105549-
}
