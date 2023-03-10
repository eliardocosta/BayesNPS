% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rpos.nps.R
\name{rpos.nps}
\alias{rpos.nps}
\title{Random generation for the posterior distribution of the NPS}
\usage{
rpos.nps(n, x, a = rep(1, 3))
}
\arguments{
\item{n}{number of observations for random generation.}

\item{x}{a vector with the observed counts of detractors, passives and promoters, respectively.}

\item{a}{vector parameter of the prior Dirichlet distribution.}
}
\value{
A vector with a random sample from the posterior distribution of the NPS.
}
\description{
Random generation for the posterior distribution of the NPS
}
\examples{
rpos.nps(n = 100, x = c(10, 50, 80))
}
