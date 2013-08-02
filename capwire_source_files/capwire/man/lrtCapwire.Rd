\name{lrtCapwire}
\alias{lrtCapwire}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood Ratio Test Comparing TIRM to ECM
}
\description{
Performs Likelihood Ratio Test to compare model fits between the Equal Capture Model (ECM) and Two Innate-Rates Model (TIRM).

A parametric bootstrap is used to generate a null-distribution for the likelihood ratio test-statistic to assess significance
}
\usage{
lrtCapwire(ecm, tirm, bootstraps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ecm}{
An object inherited from \code{fitEcm}
}
  \item{tirm}{
An object inherited from \code{fitTirm}
}
  \item{bootstraps}{
The number of parameteric bootstraps to perform to generate null-distribution for likelihood ratio
}
}
\details{
The maximum likelihood estimate for population size under the ECM model is used to simulate data under
the ECM model. The data is then fit to both the ECM and TIRM models and the likelihood ratio test is evaluated.
The observed Likelihood Ratio statistic is then compared to the distribution of likelihood ratio statistics observed when
the null model is true in order to estimate a p-value.

Note that the likelihood ratio test should not be performed on an object inherited from \code{fitTirmPartition} as
the partitioning algorithm assumes that there is heterogeneity in capture rates among individuals and that the rates
can be divided into classes

In most applications of likelihood ratio tests, the more parameterized model will always have a higher likelihood. 
However, this is not necessarily the case for the models implemented in \code{capwire}.
The likelihood involves a combinatorial term which denotes the number of ways a population
could give rise to the observed vector of counts. As such, for some data sets (such as those for which the assumptions of the ECM model are valid),
there will be more ways to obtain the observed data if the population is not subdivided into classes.
}
\value{
\item{LR}{Likelihood Ratio}
\item{p.value}{P value for the siginficance of the likelihood ratio test derived from a parametric bootstrap}

}
\references{
Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005.

Pennell M.W., C.R. Stansbury, L.P. Waits and C.R. Miller. submitted. capwire: A R Package for Estimating Population Census Size from Non-Invasive Genetic Sampling

Stansbury et al. in prep
}
\author{
Matthew W. Pennell
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitEcm}}, \code{\link{fitTirm}}, \code{\link{simEcm}}
}
\examples{
## Simulate data under Equal Capture Model

data <- simEcm(n=20, s=60)

## Fit Equal Capture Model (ECM) to count data

ecm <- fitEcm(data, max.pop=200)

## Fit Two-Innate Rates Model (TIRM) to count data

tirm <- fitTirm(data, max.pop=200)

## Perform Likelihood Ratio Test

lrtCapwire(ecm=ecm, tirm=tirm, bootstraps=10)

}

