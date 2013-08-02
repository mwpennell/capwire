\name{fitTirmPartition}
\alias{fitTirmPartition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit the Two Innate Rates Model After Partiioning Data
}
\description{
Partitions count data into three groups based on PART algorithm of Stansbury et al. in prep. The lower two groups are retained for the population estimation
procedure implemented in capwire. The upper group is excluded from the estimation procedure. This is done to handle overdispersed data which violate the assumptions of 
the Two-Innate Rates Model (TIRM).

The function also returns a p-value to evaluate whether partitioning the data is statistically warranted.

The function then fits the Two-Innate Rates Model (TIRM) to the lower set of the partitioned data and adds the individuals from the upper set to the estimate of population size
}
\usage{
fitTirmPartition(data, max.pop, n.boots.part=100)
}
\arguments{
  \item{data}{
A two-column data frame with the first column specifiying the capture class (i.e. individuals in class i were caught i times) and the second column specifying the number of individuals in each class}
  \item{max.pop}{
The maximum population size}
  \item{n.boots.part}{
The number of bootraps used to generate the null distribution for partitioning 
in order to estimate a p-value associated with the partitioning scheme used (see details).
For most applications the default value (\code{n.boots.part=100}) should be appropriate.
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

The TIRM model is then fit to the lower two classes obtained from partitioning the data in order to obtain the MLE for population size. 
Individuals from the upper class (i.e. those that were captured most often) to the population size estimate. It is therefore assumed that we have captured all individuals in the population with the highest probability of capture
and that the class designation has been made without error.

The TIRM model fit by this function assumes that individuals can be assigned to two classes. Class A represent the frequently captured individuals. Class B represents the infrequently captured individuals.

The value is specified for \code{max.pop} is not likely to matter as long as it is much greater than the maximum likelihood estimate for population size.

Note that if the data contains only singletons, the data is not informative and the maximum likelihood estimate for population size will be equal to \code{max.pop}

}
\value{
\item{model}{The model specified}
\item{likelihood}{The likelihood of the model}
\item{ml.pop.size}{The maximum likelihood estimate for population size including the addition of individuals excluded by partitioning method}
\item{ml.na}{The maximum likelihood estimate for the number of individuals in class A}
\item{ml.nb}{The maximum likelihood estimate for the number of individuals in class B}
\item{alpha}{The ratio of the rates of captures between class A and class B individuals}
\item{cap.ind}{The mean number of captures per individual}
\item{sampled.ind}{Total number of individuals in the sample (prior to partitioning)}
\item{sample.size}{Total number of samples in the data set (post partitioning)}
\item{max.pop}{The maximum population size specified by \code{max.pop}}
\item{excluded}{A two-column data frame specifying the capture classes and the number of individuals in this capture class for individuals excluded from the estimation procedure. These individuals are added to the population size estimate}
\item{p.value.partition}{The p-value obtained from the parametric bootstraps for partitioning the data into the three groups}
}

\references{
Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005.

Pennell M.W., C.R. Stansbury, L.P. Waits and C.R. Miller. submitted. capwire: A R Package for Estimating Population Census Size from Non-Invasive Genetic Sampling

Stansbury C.R., D.E. Ausband, P. Zager, C.M. Mack, C.R. Miller, M.W. Pennell, and L.P. Waits. in prep. Non-invasive genetic sampling of rendezvous sites and population estimation of grey wolves in Idaho, USA
}
\author{
Matthew W. Pennell
}

\seealso{
\code{\link{fitTirm}}, \code{\link{partitionCountData}}
}

\examples{
## Use dummy data set with a few individuals having very high capture counts
d <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,4,4,6,6,7,8,10,10,14,17,19,22,22,25)

data <- buildClassTable(d)

## Partition the data
## Fit the TIRM to the lower two classes obtained from partitioning data set

res <- fitTirmPartition(data=data, max.pop=200, n.boots.part=10)

res

## Compare population estimate from fitting TIRM 
## to partitioned data set to population estimate obtained from fitting all data

res.all <- fitTirm(data=data, max.pop=200)$ml.pop.size

res.all

}
