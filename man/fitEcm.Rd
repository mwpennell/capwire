% Generated by roxygen2 (4.0.1): do not edit by hand
\name{fitEcm}
\alias{fitEcm}
\title{Fit Equal Capture Model (ECM)}
\usage{
fitEcm(data, max.pop)
}
\arguments{
\item{data}{A two-column data frame with the first column specifiying the capture class (i.e. individuals in class i were caught i times) and the second column specifying the number of individuals in each class}

\item{max.pop}{The maximum population size}
}
\value{
an object of class \code{capfit} containing the follow elements:
\describe{
 \item{model}{The model specified}
 \item{likelihood}{The likelihood of the model}
 \item{ml.pop.size}{The maximum likelihood estimate for population size}
 \item{cap.ind}{The mean number of captures per individual}
 \item{sampled.ind}{The total number of individuals in the sample}
 \item{sample.size}{Total number of samples in the data set}
 \item{max.pop}{The maximum population size specified by \code{max.pop}}
}
}
\description{
Fits the Equal Capture Model (ECM) to count data to obtain the MLE for population size
}
\details{
The ECM model fit by this function assumes that all individuals are equally likely to be "captured".

The value is specified for \code{max.pop} is not likely to matter as long as it is much greater than the maximum likelihood estimate for population size.

Note that if the data contains only singletons, the data is not informative and the maximum likelihood estimate for population size will be equal to \code{max.pop}
}
\examples{
## Simulate data under Equal Capture Model

data <- simEcm(n=40, s=150)

## Fit Equal Capture Model to Data

res <- fitEcm(data=data, max.pop=200)

res
}
\author{
Matthew W. Pennell
}
\references{
Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005.

Pennell, M. W., C. R. Stansbury, L. P. Waits and C. R. Miller. 2013. Capwire: a R package for estimating population census size from non-invasive genetic sampling. Molecular Ecology Resources 13:154-157.
}
\seealso{
\code{\link{simEcm}}
}

