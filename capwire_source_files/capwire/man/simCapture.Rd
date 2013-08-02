\name{simCapture}
\alias{simCapture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate data under various distributions of capture rates
}
\description{
Simulates capture count data where individual capture rates are assumed to be drawn from a specified distribution.

Data can be used as input for fitting Equal Capture Model (with \code{fitEcm}), Two-Innate Rates Model (with \code{fitTirm}) or the Two-Innate Rates Model with partitioned data (with \code{fitTirmPartition})
}
\usage{
simCapture(n, s, dist.func, return.cap.probs=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
The number of individuals in the population
}
  \item{s}{
The total number of samples collected
}
  \item{dist.func}{
The distribution of capture rates within the population (see details)
}
  \item{return.cap.probs}{
Logical, signifying whether individual capture probabilities should be returned (\code{return.cap.probs=TRUE})
or not (\code{return.cap.probs=FALSE})
}
}
\details{
We assume that there is heterogeneity in the capturabilities of individuals within a population. That is, some individuals are more likely to be captured than others

We also assume that the individual capturabilities are drawn from some distribution.

The distribution is specified by the \code{dist.func} argument. \code{dist.func} takes a function with parameter n, where n specifies the number of samples to be drawn.
\code{simCapture} can take any distribution of this form but the \code{capwire} package includes several functions which allow for users to draw capture rates from several standard distributions
such as a uniform (\code{\link{drawCapRatesUnif}}), exponential (\code{\link{drawCapRatesExp}}), gamma (\code{\link{drawCapRatesGamma}}) and geometric (\code{\link{drawCapRatesGeom}}).
}
\value{
If \code{return.cap.probs=FALSE}: A two-column data frame with the first column specifiying the capture class (i.e. individuals caught i times) and the second column specifying the number of individuals in each class.

If \code{return.cap.probs=TRUE}: A list consisting of two items 1) A three-column data frame with the first column specifiying the individual IDs, the second column specifying the number of times each individual was captured
and the third denoting the capture probability associated with each individual and
2) A two-column data frame with the first column specifiying the capture class (i.e. individuals caught i times) and the second column specifying the number of individuals in each class.
}

\references{
Pennell M.W., C.R. Stansbury, L.P. Waits and C.R. Miller. submitted. capwire: A R Package for Estimating Population Census Size from Non-Invasive Genetic Sampling
}

\author{
Matthew W. Pennell
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{drawCapRatesUnif}}, \code{\link{drawCapRatesGeom}}, \code{\link{drawCapRatesExp}}, \code{\link{drawCapRatesBeta}}, \code{\link{drawCapRatesGamma}}
}

\examples{
## Specify that capture probabilities are drawn from a uniform distribution

rates <- drawCapRatesUnif(lower=0.1, upper=1)

## Simulate data using the specified distribution

data <- simCapture(n=20, s=100, dist.func=rates, return.cap.probs=FALSE)

data

## Example of sampling capture probabilities from a non-standard distribution
## We want to use a truncated normal for this purpose

## Create a function which takes argumens mean, sd, and trunc.point
drawCapRates.truncnorm <- function(mean, sd, trunc.point){
	
	## Draw a large number of samples from a normal distribution with specified mean and sd
	normal.dist <- rnorm(100000, mean, sd)
	
	## Truncate the distribution at the point specified by trunc.point
	trunc.dist <- normal.dist[which(normal.dist >= trunc.point)]
	
	## Create a function which draws n samples from trunc.dist
	function(n){
		sample(trunc.dist, size=n, replace=TRUE)
	}
	
}

## Specify the distribution

my.dist <- drawCapRates.truncnorm(mean=2, sd=1, trunc.point=0)

## Look at the function that is returned

my.dist

## Simulate data under the specified distribution

my.data <- simCapture(n=20, s=100, dist.func=my.dist, return.cap.probs=TRUE)

my.data

}
