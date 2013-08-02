\name{buildIndTable}
\alias{buildIndTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert capture counts to data frame with individuals
}
\description{
Converts a vector of capture counts into a two-column data frame consisting of all individuals and the capture counts associated with each
}
\usage{
buildIndTable(counts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
A vector of capture count data
}
}

\value{
A two-column data frame with the first column specifiying the individual IDs and the second column specifying the number of times each individual was captured.
Note here that the individual IDs are arbitrarily defined.
}
\references{
Pennell M.W., C.R. Stansbury, L.P. Waits and C.R. Miller. submitted. capwire: A R Package for Estimating Population Census Size from Non-Invasive Genetic Sampling 
}
\author{
Matthew W. Pennell
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{buildClassTable}}, \code{\link{indToClass}}, \code{\link{classToInd}}
}
\examples{
## create a vector of capture counts

mycounts <- c(1,1,1,1,1,2,2,3,3,4,5)

## build table 

data <- buildIndTable(mycounts)
data
}
