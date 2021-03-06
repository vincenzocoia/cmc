% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/YgX.R
\name{QYgX}
\alias{QYgX}
\alias{FYgX}
\title{Quantile Function for Y|X}
\usage{
QYgX(tau, ucond, cops, cpars, QY = identity)

FYgX(y, ucond, cops, cpars, FY = identity)
}
\arguments{
\item{tau}{Vector of quantile levels to evaluate at. Or,
different vectors (of the same length) can be specified
for each observation by making this a matrix (with rows corresponding
to observations).}

\item{ucond}{Matrix of uniform independent predictors. The i'th
column corresponds to the PIT score of the i'th predictor linked to the
response, conditional on previously linked predictors. See
\code{pcondseq} in the \code{copsupp} package for this.}

\item{cops}{Vector of copula model names, corresponding to the columns of
\code{ucond}, that link (predictor, response) conditional on previous
predictors in the pairing order.}

\item{cpars}{List of parameter vectors corresponding to \code{cops}.}

\item{QY, FY}{Vectorized (marginal) quantile function/cdf of the response.}

\item{y}{Vector of values to evaluate the cdf at.}
}
\value{
A matrix, with columns being the quantile levels, and rows
corresponding to the observations in \code{ucond}.
}
\description{
Computes the conditional quantile function (\code{QYgX})
or distribution function (\code{FYgX}) of a response
at specified quantile levels, given some observed predictors.
This is like the quantile version of the \code{pcondrvine} function
in the \code{copsupp} packages, except that it (1) only
takes the response to be the last entry in the vine, and (2) saves
computation time by accepting the predictors converted to an independent
set (thus bypassing the need to specify the distribution of the predictors).
}
\note{
If copulas are not permutation-symmetric, ensure that the copula
specified has its fist argument corresponding to the predictor, and second
argument the response.
}
\examples{
(dat <- matrix(runif(3*6), ncol = 3))
tau <- c(0.1, 0.3, 0.7, 0.9)
cops <- c("gum", "bvtcop", "frk")
cpars <- list(3.5, c(0.6, 3), 2.5)

## Get tau quantiles for the predictors in `dat`
QYgX(tau, dat, cops = cops, cpars = cpars, QY = qexp)

## Get sample for these predictors:
(v <- runif(nrow(dat)))
(y <- QYgX(matrix(v), dat, cops = cops, cpars = cpars, QY = qexp))
FYgX(y, dat, cops, cpars, FY = pexp)  # Same as v.
}
