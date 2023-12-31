% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betpp.R
\name{betpp}
\alias{betpp}
\title{Posterior predictive based Bayesian multigroup equivalence testing}
\usage{
betpp(values, groups, em, A, B = 10000)
}
\arguments{
\item{values}{A vector of measurements sorted in the same order as the \code{groups} variable.}

\item{groups}{A vector of groups labels corresponding to the individual measurements in the \code{groups} variable.}

\item{em}{A c x 2 matrix of lower and upper equivalence margins. Here, c is the number of pairwise comparisons of interest.}

\item{A}{A c x k matrix of pairwise contrasts. Here, k is the number of groups, i.e., \code{length(unique(groups))}.}

\item{B}{A positive integer specifying the number of posterior predictive samples to draw. By default \code{B} is set to 10000.}
}
\value{
The function returns a list object containing the following:
\itemize{
\item prob: The probability that future differences fall within the equivalence margins.
\item delta: A B x c matrix of posterior predictive samples of future differences for each pairwise comparison of interest.
}
}
\description{
Function provides the necessary tools to carry out Bayesian multigroup equivalence testing based on sampling of the posterior predictive distribution.
The function returns posterior predictive samples of future differences amongst groups.
}
\examples{
### Multigroup equivalence test for A vs. B and A vs. C
values <- rnorm(75)
groups <- rep(LETTERS[1:3], each = 25)

mad1 <- 0.65  # The equivalence margin for A vs. B
mad2 <- 0.65  # The equivalence margin for A vs. C
mads <- c(mad1, mad2)
mads <- cbind(-mads, mads)

A <- apc(3)
A <- A[1:2, ]

out <- betpp(values, groups, mads, A, B = 10000)

out$prob   # The probability that future differences
           # fall within the equivalence margins

}
\references{
Pourmohamad, T. and Lee, H.K.H. (2023). Equivalence Testing for Multiple Groups. Stat, e645.
}
