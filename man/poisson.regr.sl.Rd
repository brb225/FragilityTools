% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{poisson.regr.sl}
\alias{poisson.regr.sl}
\title{Calculate level of smallest HDR which contains the modified observation}
\usage{
poisson.regr.sl(y, ymod, hatbeta, Xregr)
}
\arguments{
\item{y}{a numeric vector of outcomes}

\item{hatbeta}{a numeric vector of estimated coefficients}

\item{Xregr}{a numeric matrix of covariates used to explain the outcomes}

\item{y.mod}{a numeric vector of outcomes after the fragility index modification}
}
\value{
a scalar numeric: a ratio of the conditional likilihoods in y and y.mod, estimated
using the MLE with the data y
}
\description{
Used for internal functions in the study of generalized fragility indices.
}
\examples{
y <- rpois(100, 1)
y.mod <- rpois(100, 1.2) # should be returned from an FI procedure
FragilityTools:::poisson.regr.sl(y, y.mod, 0, matrix(rnorm(100), ncol=100))
}
