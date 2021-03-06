% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/from_cnqr-dat2udat.R
\name{dat2udat}
\alias{dat2udat}
\title{Convert raw data to uniform}
\usage{
dat2udat(dat, cdf)
}
\arguments{
\item{dat}{A data frame or matrix of raw data.}

\item{cdf}{A list of vectorized distribution functions, corresponding
to the columns of dat.}
}
\value{
A matrix of uniform scores with the same dimensions as \code{dat},
obtained by applying the cdf functions to the data.
}
\description{
Converts a data frame/matrix into uniform PIT scores.
}
