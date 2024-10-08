% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ParamTransfoCompRisks.R
\name{likIFG.cmprsk.Cholesky}
\alias{likIFG.cmprsk.Cholesky}
\title{Full likelihood (including estimation of control function).}
\usage{
likIFG.cmprsk.Cholesky(
  parhatG,
  data,
  eoi.indicator.names,
  admin,
  conf,
  Zbin,
  inst
)
}
\arguments{
\item{parhatG}{The full parameter vector.}

\item{data}{Data frame.}

\item{eoi.indicator.names}{Vector of names of the censoring indicator columns
pertaining to events of interest. Events of interest will be modeled allowing
dependence between them, whereas all censoring events (corresponding to
indicator columns not listed in \code{eoi.indicator.names}) will be treated
as independent of every other event. If \code{eoi.indicator.names == NULL},
all events will be modelled dependently.}

\item{admin}{Boolean value indicating whether the data contains
administrative censoring.}

\item{conf}{Boolean value indicating whether the data contains confounding
and hence indicating the presence of Z and W.}

\item{Zbin}{Boolean value indicating whether the confounding variable is
binary.}

\item{inst}{Type of instrumental function to be used.}
}
\value{
Full model log-likelihood evaluation.
}
\description{
This function defines the 'full' likelihood of the model.
Specifically, it includes the estimation of the control function in the
computation of the likelihood. This function is used in the estimation of the
variance of the estimates (variance.cmprsk.R).
}
