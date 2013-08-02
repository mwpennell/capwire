\name{drawCapRatesGamma}
\alias{drawCapRatesGamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw Capture Rates from Gamma Distribution
}
\description{
This function is just a wrapper for \code{rgamma} which can be used in conjuction with \code{simCapture} to simulate capture data where the rates are drawn from a gamma distribution
}
\usage{
drawCapRatesGamma(shape, rate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{shape}{
Shape parameter of gamma distribution
}
  \item{rate}{
Rate parameter of gamma distribution
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{simCapture}}, \code{\link{drawCapRatesUnif}}, \code{\link{drawCapRatesExp}}, \code{\link{drawCapRatesBeta}}, \code{\link{drawCapRatesGeom}}
}

\examples{
## Specify the distribution

dist <- drawCapRatesGamma(shape=0.5, rate=0.5)

## Simulate a data set with the capture rates drawn from dist

data <- simCapture(n=30, s=100, dist.func=dist)

data
}
