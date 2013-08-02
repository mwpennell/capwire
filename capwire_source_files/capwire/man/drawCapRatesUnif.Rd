\name{drawCapRatesUnif}
\alias{drawCapRatesUnif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw Capture Rates from Uniform Distribution
}
\description{
This function is just a wrapper for \code{runif} which can be used in conjuction with \code{simCapture} to simulate capture data where the rates are drawn from a uniform distribution
}
\usage{
drawCapRatesUnif(lower, upper)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lower}{
Lower bound of the unifrom distribution
}
  \item{upper}{
Upper bound of the uniform distribution
}
}
\details{
This function returns a function which can be used as an argument for \code{simCapture}
}
\value{
\item{function}{A function which takes the number of capture rates to draw as an argument}
}
\references{
Pennell M.W., C.R. Stansbury, L.P. Waits and C.R. Miller. submitted. capwire: A R Package for Estimating Population Census Size from Non-Invasive Genetic Sampling
}
\author{
Matthew W. Pennell
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{simCapture}}, \code{\link{drawCapRatesGeom}}, \code{\link{drawCapRatesExp}}, \code{\link{drawCapRatesBeta}}, \code{\link{drawCapRatesGamma}}
}
\examples{
## Specify the distribution

dist <- drawCapRatesUnif(lower=0.01, upper=1)

## Simulate a data set with the capture rates drawn from dist

data <- simCapture(n=30, s=100, dist.func=dist)

data
}