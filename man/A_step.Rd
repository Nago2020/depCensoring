% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoundingCovariateEffects.R
\name{A_step}
\alias{A_step}
\title{A-step in the EAM algorithm described in KMS19}
\usage{
A_step(evaluations, verbose = 0)
}
\arguments{
\item{evaluations}{Matrix containing each point that was already evaluated,
alongside the corresponding test statistic and critical value, as its rows.}

\item{verbose}{Verbosity parameter.}
}
\value{
Results of the A-step.
}
\description{
This function performs the approximation step in the EAM
algorithm. More specifically, it fits a Gaussian-process regression model
(Kriging) to the evaluated data points \eqn{(\theta, c(\theta))}.
}
\seealso{
Package \pkg{rkriging}.
}
