\name{drawCapRatesExp}
\alias{drawCapRatesExp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw Capture Rates from Exponential Distribution
}
\description{
This function is just a wrapper for \code{rexp} which can be used in conjuction with \code{simCapture} to simulate capture data where the rates are drawn from an exponential distribution
}
\usage{
drawCapRatesExp(r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
Rate parameter of an exponenital distirbution
}
}
\details{
This function returns a function which can be used as an argument for \code{simCapture}
}
\value{
A function which takes the number of capture rates to draw as an argument
}
\references{
Pennell M.W., C.R. Stansbury, L.P. Waits and C.R. Miller. submitted. capwire: A R Package for Estimating Population Census Size from Non-Invasive Genetic Sampling
}
\author{
Matthew W. Pennell
}


\seealso{
\code{\link{simCapture}}, \code{\link{drawCapRatesUnif}}, \code{\link{drawCapRatesGamma}}, \code{\link{drawCapRatesBeta}}, \code{\link{drawCapRatesGeom}}
}

\examples{
## Specify the distribution

dist <- drawCapRatesExp(r=0.5)

## Simulate a data set with the capture rates drawn from dist

data <- simCapture(n=30, s=100, dist.func=dist)

data
}
