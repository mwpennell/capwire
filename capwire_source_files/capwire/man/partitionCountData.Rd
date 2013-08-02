\name{partitionCountData}
\alias{partitionCountData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partition the Count Data
}
\description{
Partitions count data into three groups based on PART algorithm of Stansbury et al. in prep. The lower two groups are retained for the population estimation
procedure implemented in capwire. The upper group is excluded from the estimation procedure. This is done to handle overdispersed data which violate the assumptions of 
the Two-Innate Rates Model (TIRM).

The function also returns a p-value to evaluate whether partitioning the data is statistically warranted.
}
\usage{
partitionCountData(data, n.boots.part = 100, max.pop=500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A two-column data frame with the first column specifying the individual IDs and the second column specifying the number of times each individual was identified
}
  \item{n.boots.part}{
Number of parametric bootstraps to perform to obtain a null distribution of the test-statisic for partitioning the data (see details)
}
  \item{max.pop}{
The maximum population size used for fitting the Two-Innate Rates Model (see \code{\link{fitTirm}} for details)
}  
}
\details{
This function implements a new partitioning (PART) method that removes individuals detected a large number of times from the data. 
These individuals provide little information about the population size and, because they are inconsistent with the modeling assumption of just two capture rates, they cause estimation problems in capwire.
Specifically, the PART algorithm attempts to remove very large capture counts and leave individuals captures between few and a moderate number of times.

The PART algorithm first fits the Two-Innate Rates Model (TIRM) to the entire data set.
Assuming the distribution of capture counts can be approximated as either a two (null) or three (alternative) class multinomial distribution, the
algorithm considers all possible ways to partition the count data under both the null and alternative distributions. The partitioning schemes which
maximizes the respective likelihoods is retained. A test-statistic, the ratio of multinomial (log) likelihoods (3:2) is obtained for the observed data.

The MLE for number of individuals in class a, number of individuals in class b and the ratio of capture probabilities (alpha) are used to simulate \code{n.boots.part} data sets
under the TIRM. The test-statistic is calculated for each data set as described above. 

The observed test-statistic is then compared to the distribution of simulated test-statistics. When the observed test-statistic falls in the tail of the distribution of the test-statisitic from the simulated data, the null model can be rejected.
}
\value{
\item{low.2.third}{A two-column data frame specifying the capture classes and the number of individuals in this capture class for the lower group}
\item{up.1.third}{A two-column data frame specifying the capture classes and the number of individuals in this capture class for the upper group. This group will be excluded from the estimation procedure in \code{\link{fitTirmPartition}}}
\item{p.2.3}{The p-value obtained from the parametric bootstraps for partitioning the data into the three groups}
}
\references{
Pennell M.W., C.R. Stansbury, L.P. Waits and C.R. Miller. submitted. capwire: A R Package for Estimating Population Census Size from Non-Invasive Genetic Sampling

Stansbury C.R., D.E. Ausband, P. Zager, C.M. Mack, C.R. Miller, M.W. Pennell, and L.P. Waits. in prep. Non-invasive genetic sampling of rendezvous sites and population estimation of grey wolves in Idaho, USA
}
\author{
Craig R. Miller and Matthew W. Pennell
}
\note{
This method is described in detail in Stansbury et al. in prep

The function \code{partitionCountData} is used internally in \code{\link{fitTirmPartition}}
}


\seealso{
\code{\link{fitTirmPartition}}, \code{\link{fitTirm}}
}
\examples{
## Use dummy data set with a few individuals having very high capture counts
d <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,4,4,6,6,7,8,10,10,14,17,19,22,22,25)

data <- buildClassTable(d)

## Perform Partitioning Algorithm

part.data <- partitionCountData(data=data, n.boots.part=10, max.pop=500)

part.data

}
