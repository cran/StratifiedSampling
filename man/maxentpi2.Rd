% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{maxentpi2}
\alias{maxentpi2}
\title{Joint inclusion probabilities of maximum entropy.}
\usage{
maxentpi2(pikr)
}
\arguments{
\item{pikr}{A vector of inclusion probabilities.}
}
\value{
A matrix, the joint inclusion probabilities.
}
\description{
This function computes the matrix of the joint inclusion of the maximum entropy sampling with fixed sample size. It can handle unequal inclusion probabilities.
}
\details{
The sampling design maximizes the entropy design:
\deqn{I(p) = - \sum s p(s) log[p(s)].}

This function is a C++ implementation of \code{\link[sampling:UPMEpik2frompikw]{UPMEpik2frompikw}}.
More details could be find in Tille (2006).
}
\references{
Tille, Y. (2006), Sampling Algorithms, springer
}
\author{
Raphaël Jauslin \email{raphael.jauslin@unine.ch}
}
