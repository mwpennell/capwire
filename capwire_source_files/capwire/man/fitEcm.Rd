\name{fitEcm}
\alias{fitEcm}
\title{
Fit Equal Capture Model (ECM)
}
\description{
Fits the Equal Capture Model (ECM) to count data to obtain the MLE for population size
}
\usage{
fitEcm(data, max.pop)
}
\arguments{
  \item{data}{
A two-column data frame with the first column specifiying the capture class (i.e. individuals in class i were caught i times) and the second column specifying the number of individuals in each class}
  \item{max.pop}{
The maximum population size}
}
\details{
The ECM model fit by this function assumes that all individuals are equally likely to be "captured". 

The value is specified for \code{max.pop} is not likely to matter as long as it is much greater than the maximum likelihood estimate for population size.

Note that if the data contains only singletons, the data is not informative and the maximum likelihood estimate for population size will be equal to \code{max.pop}
}
\value{
\item{model}{The model specified}
\item{likelihood}{The likelihood of the model}
\item{ml.pop.size}{The maximum likelihood estimate for population size}
\item{cap.ind}{The mean number of captures per individual}
\item{sampled.ind}{The total number of individuals in the sample}
\item{sample.size}{Total number of samples in the data set}
\item{max.pop}{The maximum population size specified by \code{max.pop}}
}
\references{
Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005.
}
\author{
Matthew W. Pennell
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{simEcm}}
}
\examples{
## Simulate data under Equal Capture Model

data <- simEcm(n=40, s=150)

## Fit Equal Capture Model to Data

res <- fitEcm(data=data, max.pop=200)

res
}
