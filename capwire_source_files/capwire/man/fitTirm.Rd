\name{fitTirm}
\alias{fitTirm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Two Innate Rates Model (TIRM)
}
\description{
Fits the Two Innate Rates Model (TIRM) to count data to obtain the MLE for population size
}
\usage{
fitTirm(data, max.pop)
}
\arguments{
  \item{data}{
A two-column data frame with the first column specifiying the capture class (i.e. individuals in class i were caught i times) and the second column specifying the number of individuals in each class}
  \item{max.pop}{
The maximum population size}
}
\details{
The TIRM model fit by this function assumes that individuals can be assigned to two classes. Class A represent the frequently captured individuals. Class B represents the infrequently captured individuals.

The value is specified for \code{max.pop} is not likely to matter as long as it is much greater than the maximum likelihood estimate for population size.

Note that if the data contains only singletons, the data is not informative and the maximum likelihood estimate for population size will be equal to \code{max.pop}
}
\value{
\item{model}{The model specified}
\item{likelihood}{The likelihood of the model}
\item{ml.pop.size}{The maximum likelihood estimate for population size}
\item{ml.na}{The maximum likelihood estimate for the number of individuals in class A}
\item{ml.nb}{The maximum likelihood estimate for the number of individuals in class B}
\item{alpha}{The ratio of the rates of captures between class A and class B individuals}
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
\code{\link{simTirm}}
}

\examples{
## Simulate data under Two Innate Rates Model

data <- simTirm(na=20, nb=15, alpha=4, s=150)

## Fit Two Innate Rates Model to Data

res <- fitTirm(data=data, max.pop=200)

res
}
