% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sibcorr.R
\name{sibcorr}
\alias{sibcorr}
\title{Estimate sibling and cousin correlations.}
\usage{
sibcorr(formula, data, weight = 4, restriction = NULL,
  variance = "pre", reps = 0, ci_level = 95, cores = 1)
}
\arguments{
\item{formula}{Three-part formula describing the outcome variable, control variables to be regressed out, and identifiers for the individual, immediate family, and extended family}

\item{data}{Estimation data set.}

\item{weight}{Select one of four weighting schemes.}

\item{restriction}{Put a restriction on the pairs to be used for estimating the covariance}

\item{variance}{Estimate variance on all individuals ("pre", the default) or only restricted sample ("post")}

\item{reps}{Number of bootstrap replications - if set to zero, confidence intervals are not estimated}

\item{ci_level}{Set level for bootstrap confidence interval}

\item{cores}{Set number of processor cores to use for bootstrap estimation}
}
\value{
The estimated sibling or cousin correlation coefficient,
number of individuals, sibling or cousin pairs, and families or extended families,
and optionally a bootstrap confidence interval
}
\description{
\code{sibcorr} estimates sibling or cousin correlations with optional block bootstrap standard errors.
}
\details{
The formula must be specified as \code{outcome ~ controls | individual + family + ext_family},
where \code{individual} is an individual identifier, \code{family} is a family (sibling group) identifier,
and \code{ext_family} is an extended family (cousin group) identifier.
If \code{ext_family} is omitted, the sibling correlation is estimated -
otherwise, the cousin correlation is estimated

The formula does not handle functions on the left hand side.
This means that any transformations of the outcome variable must be performed before estimation.

The restriction is specified as a vector with the first element giving the name of the variable to restrict on.
If the second element is 0, only pairs with the same value are used.
If the second element is a positive integer, only pairs with that specific difference are used.
If the second element is "unequal", only pairs with different values for the variable are used.
If two integers are given (as second and third elements), and the first is smaller than the second,
only pairs with a difference within that range (inclusive) are used;
if the second is smaller than the first,
only pairs with a difference outside that range (exclusive) are used.
}
