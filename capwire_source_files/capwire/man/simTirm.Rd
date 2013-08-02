\name{simTirm}
\alias{simTirm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Data Under Two-Innate Rates Model
}
\description{
Simulates capture count data under the assumptions of the Two-Innate Rates Model (TIRM) where the individuals in the population are assumed to come from two classes: easy to capture individuals and difficult to capture individuals.
}
\usage{
simTirm(na, nb, alpha, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{na}{
Number of individuals in the "easy to capture" class  
}
  \item{nb}{
Number of individuals in the "difficult to capture" class
}
  \item{alpha}{
The ratio of the capture rates between the two classes (capture rate of class a individual/capture rate of class b individuals)
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
\code{\link{fitTirm}}
}
\examples{
## Simualte data under Two-Innate Rates Model

data <- simTirm(na=20, nb=10, alpha=4, s=100)

## Fit the Two-Innate Rates Model 

tirm <- fitTirm(data=data, max.pop=200)

tirm
}