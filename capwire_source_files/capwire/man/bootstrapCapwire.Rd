\name{bootstrapCapwire}
\alias{bootstrapCapwire}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parametric Bootstrap ECM and TIRM Models
}
\description{
Uses maximum likelihood parameter estimates from \code{fitEcm}, \code{fitTirm}, or \code{fitTirmPartition} to perform a parametric bootstrap to get confidence intervals for the estimate of the population size
}
\usage{
bootstrapCapwire(x, bootstraps = 1000, CI = c(0.025, 0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object inherited from \code{\link{fitEcm}}, \code{\link{fitTirm}}, or \code{\link{fitTirmPartition}}
}
  \item{bootstraps}{
The number of bootstraps to be performed (default is 1000)
}
  \item{CI}{
A vector of quantiles to to generate a confidence interval for the population estimate. The default is \code{c(0.025, 0.975)}, 
denoting a symmetrical 95 percent confidence interval.
}
}
\details{
This function uses the ML estimates obtained from fitting the model to simulate data under the model.

\code{\link{bootstrapCapwire}} inherits an object from \code{\link{fitEcm}}, \code{\link{fitTirm}}, or \code{\link{fitTirmPartition}} such that the model and parameter estimates do not need to be specified.

The ML estimate for the population size will also be returned but this will not be changed by \code{bootstrapCapwire}

Note that if the model is a poor fit to the data, the confidence intervals may not be reliable.

The lower confidence interval is bounded by the number of unique individuals in the sample.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{ml.pop.size}{The maximum likelihood estimate for the population size obtained by fitting the model}
\item{conf.int}{The confidence interval for the estimate of the population size}
%% ...
}
\references{
Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005.
}
\author{
Matthew W. Pennell
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitEcm}}, \code{\link{fitTirm}}, \code{\link{fitTirmPartition}}
}
\examples{
## Simulate data under Equal Capture Model

data <- simEcm(n=40, s=150)

## Fit Equal Capture Model to Data

res <- fitEcm(data=data, max.pop=200)

## Perform bootstrap to get confidence intervals

ci <- bootstrapCapwire(x=res, bootstraps=50, CI = c(0.025, 0.975))

}
