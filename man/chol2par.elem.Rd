% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ParamTransfoCompRisks.R
\name{chol2par.elem}
\alias{chol2par.elem}
\title{Transform Cholesky decomposition to covariance matrix parameter element.}
\usage{
chol2par.elem(a, b, par.chol1)
}
\arguments{
\item{a}{The row index of the covariance matrix element to be computed.}

\item{b}{The column index of the covariance matrix element to be computed.}

\item{par.chol1}{The vector of Cholesky parameters.}
}
\value{
Specified element of the covariance matrix resulting from the
provided Cholesky decomposition.
}
\description{
This function transforms the parameters of the Cholesky de-
composition to a covariance matrix element. This function is used in
chol2par.R.
}
