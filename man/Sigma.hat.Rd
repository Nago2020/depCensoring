% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoundingCovariateEffects.R
\name{Sigma.hat}
\alias{Sigma.hat}
\title{Compute the variance-covariance matrix of the moment functions.}
\usage{
Sigma.hat(data, beta, t, hp, m.avg = NULL, mi.mat = NULL)
}
\arguments{
\item{data}{Data frame.}

\item{beta}{Coefficient vector.}

\item{t}{Time point of interest.}

\item{hp}{List of hyperparameters.}

\item{m.avg}{A precomputed vector of the sample average of the moment
functions. If not supplied, this vector is computed. Default is
\code{m.avg = NULL}.}

\item{mi.mat}{A precomputed matrix of moment function evaluations at each
observation. If supplied, some computations can be skipped. Default is
\code{mi.mat = NULL}.}
}
\description{
This function comptutes the empricical variance-covariance
matrix of the moment functions.
}
