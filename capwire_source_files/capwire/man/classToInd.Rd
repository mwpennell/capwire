\name{classToInd}
\alias{classToInd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a table of capture class data into a table of individuals
}
\description{
Converts a data frame where individuals are grouped by capture class into a data frame consisting of individuals and their associated capture counts. 
}
\usage{
classToInd(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A two-column data frame with the first column specifiying the capture class (where all individuals in class i were caught i times) and the second column specifying the number of individuals in this capture class.
}
}
\details{
This is primarily used as an internal function through \code{capwire}
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
\code{\link{buildIndTable}}, \code{\link{buildClassTable}}, \code{\link{indToClass}}
}
\examples{
## create a vector of capture counts

counts <- c(1,1,1,2,2,2,3,3,4)

## make into a table of classes

class.data <- buildClassTable(counts)

## convert to table of individuals

data <- classToInd(class.data)

data
}
