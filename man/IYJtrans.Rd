% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yeo-johnsonTransformation.R
\name{IYJtrans}
\alias{IYJtrans}
\title{Inverse Yeo-Johnson transformation function}
\usage{
IYJtrans(y, theta)
}
\arguments{
\item{y}{The argument to be supplied to the inverse Yeo-Johnson transformation.}

\item{theta}{The parameter of the inverted Yeo-Johnson transformation. This
should be a number in the range [0,2].}
}
\value{
The transformed value of y.
}
\description{
Computes the inverse Yeo-Johnson transformation of the provided
argument.
}
