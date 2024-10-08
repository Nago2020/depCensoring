% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ParamTransfoCompRisks.R
\name{estimate.cmprsk}
\alias{estimate.cmprsk}
\title{Estimate the competing risks model of Rutten, Willems et al. (20XX).}
\usage{
estimate.cmprsk(
  data,
  admin,
  conf,
  eoi.indicator.names = NULL,
  Zbin = NULL,
  inst = "cf",
  realV = NULL,
  eps = 0.001
)
}
\arguments{
\item{data}{A data frame, adhering to the following formatting rules:
\itemize{
   \item The first column, named \code{"y"}, contains the observed times.
   \item The next columns, named \code{"delta1"}, \code{delta2}, etc. contain
   the indicators for each of the competing risks.
   \item The next column, named \code{da}, contains the censoring indicator
   (independent censoring).
   \item The next column should be a column of all ones (representing the
   intercept), names \code{x0}.
   \item The subsequent columns should contain the values of the covariates,
   named \code{x1}, \code{x2}, etc.
   \item When applicable, the next column should contain the values of the
   endogenous variable. This column should be named \code{z}.
   \item When \code{z} is provided and an instrument for \code{z} is
   available, the next column, named \code{w}, should contain the values for
   the instrument.
}}

\item{admin}{Boolean value indicating whether the data contains
administrative censoring.}

\item{conf}{Boolean value indicating whether the data contains confounding
and hence indicating the presence of \code{z} and, possibly, \code{w}.}

\item{eoi.indicator.names}{Vector of names of the censoring indicator columns
pertaining to events of interest. Events of interest will be modeled allowing
dependence between them, whereas all censoring events (corresponding to
indicator columns not listed in \code{eoi.indicator.names}) will be treated
as independent of every other event. If \code{eoi.indicator.names == NULL},
all events will be modeled dependently.}

\item{Zbin}{Indicator value indicating whether (\code{Zbin = TRUE}) or not
\code{Zbin = FALSE} the endogenous covariate is binary. Default is
\code{Zbin = NULL}, corresponding to the case when \code{conf == FALSE}.}

\item{inst}{Variable encoding which approach should be used for dealing with
the confounding. \code{inst = "cf"} indicates that the control function
approach should be used. \code{inst = "W"} indicates that the instrumental
variable should be used 'as is'. \code{inst = "None"} indicates that Z will
be treated as an exogenous covariate. Finally, when \code{inst = "oracle"},
this function will access the argument \code{realV} and use it as the values
for the control function. Default is \code{inst = "cf"}.}

\item{realV}{Vector of numerics with length equal to the number of rows in
\code{data}. Used to provide the true values of the instrumental function
to the estimation procedure.}

\item{eps}{Value that will be added to the diagonal of the covariance matrix
during estimation in order to ensure strictly positive variances.}
}
\value{
A list of parameter estimates in the second stage of the estimation
algorithm (hence omitting the estimates for the control function), as well
as an estimate of their variance and confidence intervals.
}
\description{
This function estimates the parameters in the competing risks
model described in Willems et al. (2024+). Note that this model
extends the model of Crommen, Beyhum and Van Keilegom (2024) and as such, this
function also implements their methodology.
}
\examples{
\donttest{

n <- 200

# Set parameters
gamma <- c(1, 2, 1.5, -1)
theta <- c(0.5, 1.5)
eta1 <- c(1, -1, 2, -1.5, 0.5)
eta2 <- c(0.5, 1, 1, 3, 0)

# Generate exogenous covariates
x0 <- rep(1, n)
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.5)

# Generate confounder and instrument
w <- rnorm(n)
V <- rnorm(n, 0, 2)
z <- cbind(x0, x1, x2, w) \%*\% gamma + V
realV <- z - (cbind(x0, x1, x2, w) \%*\% gamma)

# Generate event times
err <- MASS::mvrnorm(n, mu = c(0, 0), Sigma =
matrix(c(3, 1, 1, 2), nrow = 2, byrow = TRUE))
bn <- cbind(x0, x1, x2, z, realV) \%*\% cbind(eta1, eta2) + err
Lambda_T1 <- bn[,1]; Lambda_T2 <- bn[,2]
x.ind = (Lambda_T1>0)
y.ind <- (Lambda_T2>0)
T1 <- rep(0,length(Lambda_T1))
T2 <- rep(0,length(Lambda_T2))
T1[x.ind] = ((theta[1]*Lambda_T1[x.ind]+1)^(1/theta[1])-1)
T1[!x.ind] = 1-(1-(2-theta[1])*Lambda_T1[!x.ind])^(1/(2-theta[1]))
T2[y.ind] = ((theta[2]*Lambda_T2[y.ind]+1)^(1/theta[2])-1)
T2[!y.ind] = 1-(1-(2-theta[2])*Lambda_T2[!y.ind])^(1/(2-theta[2]))
# Generate adminstrative censoring time
C <- runif(n, 0, 40)

# Create observed data set
y <- pmin(T1, T2, C)
delta1 <- as.numeric(T1 == y)
delta2 <- as.numeric(T2 == y)
da <- as.numeric(C == y)
data <- data.frame(cbind(y, delta1, delta2, da, x0, x1, x2, z, w))
colnames(data) <- c("y", "delta1", "delta2", "da", "x0", "x1", "x2", "z", "w")

# Estimate the model
admin <- TRUE                # There is administrative censoring in the data.
conf <- TRUE                 # There is confounding in the data (z)
eoi.indicator.names <- NULL  # We will not impose that T1 and T2 are independent
Zbin <- FALSE                # The confounding variable z is not binary
inst <- "cf"                 # Use the control function approach
# Since we don't use the oracle estimator, this argument is ignored anyway
realV <- NULL
estimate.cmprsk(data, admin, conf, eoi.indicator.names, Zbin, inst, realV)

}

}
\references{
Willems et al. (2024+). Flexible control function approach under competing risks (in preparation).

Crommen, G., Beyhum, J., and Van Keilegom, I. (2024). An instrumental variable approach under dependent censoring. Test, 33(2), 473-495.
}
