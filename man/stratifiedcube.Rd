% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stratifiedcube.R
\name{stratifiedcube}
\alias{stratifiedcube}
\title{Stratified Sampling}
\usage{
stratifiedcube(X, strata, pik, EPS = 1e-07)
}
\arguments{
\item{X}{A matrix of size (\eqn{N} x \eqn{p}) of auxiliary variables on which the sample must be balanced.}

\item{strata}{A vector of integers that specifies the stratification..}

\item{pik}{A vector of inclusion probabilities.}

\item{EPS}{epsilon value}
}
\value{
A vector with elements equal to 0 or 1. The value 1 indicates that the unit is selected while the value 0 is for rejected units.
}
\description{
This function implements a method for selecting a stratified sample. It really improves the performance of the function \code{\link{fbs}} and \code{\link{balstrat}}.
}
\details{
The function is selecting a balanced sample very quickly even if the sum of inclusion probabilities within strata are non-integer. The function should be used in preference. Firstly, a flight phase is performed on each strata. Secondly, the function \code{\link{findB}} is used to find a particular matrix to apply a flight phase by using the cube method proposed by Chauvet, G. and Tillé, Y. (2006). Finally, a landing phase is applied by suppression of variables.
}
\examples{
N <- 100
n <- 10
p <- 4
X <- matrix(rgamma(N*p,4,25),ncol = p)
strata <- as.matrix(rep(1:n,each = N/n))
pik <- rep(n/N,N)

s <- stratifiedcube(X,strata,pik)

t(X/pik)\%*\%s
t(X/pik)\%*\%pik

Xcat <- disj(strata)

t(Xcat)\%*\%s
t(Xcat)\%*\%pik

}
\references{
Chauvet, G. and Tillé, Y. (2006). A fast algorithm of balanced sampling. \emph{Computational Statistics}, 21/1:53-62
}
\seealso{
\code{\link{fbs}}, \code{\link{balstrat}}, \code{\link{landingRM}}, \code{\link{ffphase}}
}
