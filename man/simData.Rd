% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simData.R
\name{simData}
\alias{simData}
\title{Simulate data}
\usage{
simData(
  n = 200,
  p = 10,
  L = 3,
  Sigma = 0,
  kappas = 2,
  proportion.model = "dirichlet",
  phi = 10
)
}
\arguments{
\item{n}{TBA}

\item{p}{TBA}

\item{L}{TBA}

\item{Sigma}{TBA}

\item{kappas}{value of the Weibull's shape parameter}

\item{proportion.model}{One of \code{c("alr", "cloglog", "log", "dirichlet")}}

\item{phi}{TBA}
}
\value{
An object of ...
}
\description{
Simulate survival data. (This requires external function 'metropolis_sampler')
}
\examples{

x <- 1

}
