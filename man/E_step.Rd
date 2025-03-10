% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoundingCovariateEffects.R
\name{E_step}
\alias{E_step}
\title{E-step in the EAM algorithm as described in KMS19.}
\usage{
E_step(thetas, test.fun, dir, evaluations, verbose)
}
\arguments{
\item{thetas}{Points at which to perform the E-step. Usually the result of
the M-step.}

\item{test.fun}{Function returning the test statistic, as well as the critical
value.}

\item{dir}{Direction in which to optimize. For finding upper bounds, set
\code{dir = 1}, for finding lower bounds, set \code{dir = -1}.}

\item{evaluations}{Matrix containing each point that was already evaluated,
alongside the corresponding test statistic and critical value, as its rows.}

\item{verbose}{Verbosity parameter.}
}
\value{
Results of the E-step.
}
\description{
This function performs the estimation step in the EAM algorithm.
}
