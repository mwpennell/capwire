\name{buildClassTable}
\alias{buildClassTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert capture counts to data frame with capture classes
}
\description{
Converts a vector of capture counts into a two-column data frame consisting of all capture classes and the individuals associated with each class.
}
\usage{
buildClassTable(counts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
A vector of capture count data
}
}

\value{
A two-column data frame with the first column specifiying the capture class (where all individuals in class i were caught i times) and the second column specifying the number of individuals in this capture class.

The data frame can be used as the data argument for any of the model-fitting functions implemented in capwire
}
\references{
Pennell M.W., C.R. Stansbury, L.P. Waits and C.R. Miller. submitted. capwire: A R Package for Estimating Population Census Size from Non-Invasive Genetic Sampling 
}
\author{
Matthew W. Pennell
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{buildIndTable}}, \code{\link{indToClass}}, \code{\link{classToInd}}
}
\examples{
## create a vector of capture counts

counts <- c(1,1,1,1,1,2,2,3,3,4,5)

## build table 

data <- buildClassTable(counts)
data
}
