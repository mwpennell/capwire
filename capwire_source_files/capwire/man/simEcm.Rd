\name{simEcm}
\alias{simEcm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Data Under Equal Capture Model
}
\description{
Simulates capture count data under the assumptions of the Equal Capture Model (ECM), where all individuals are assumed to have an equal probability of being captured
}
\usage{
simEcm(n, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of individuals in the population
}
  \item{s}{
Total number of samples collected
}
}

\value{
A two-column data frame with the first column specifiying the capture class (i.e. individuals caught i times) and the second column specifying the number of individuals in each class.

The data frame can be used as the data argument for any of the model-fitting functions implemented in capwire
}
\references{
Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005.

Pennell M.W., C.R. Stansbury, L.P. Waits and C.R. Miller. submitted. capwire: A R Package for Estimating Population Census Size from Non-Invasive Genetic Sampling 
}
\author{
Matthew W. Pennell
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitEcm}}
}
\examples{
## Simualte data under the Equal Capture Model

data <- simEcm(n=25, s=100)

## Fit the Equal Capture Model 

ecm <- fitEcm(data=data, max.pop=200)

ecm
}
