\name{wombat}
\alias{wombat}
\docType{data}
\title{
Wombat Capture Count Data
}
\description{
This is a data set consisting of capture count data published in Banks et al. 2003 from northern hairy-nosed wombats (Lasiorhinus kretii).
It is used as an example in the \code{capwire} tutorial.
}
\usage{
data(wombat)
}
\format{
A two-column data frame

capture.class: class assignment where individuals in class i were caught i times

No.Ind: the number of individuals in each capture.class
}
\details{
This data set can be used as input for \code{\link{fitEcm}}, \code{\link{fitTirm}}, \code{\link{fitTirmPartition}}, \code{\link{partitionCountData}}.
}
\source{
Banks, S. C., S. D. Hoyle, A. Horsup, P. Sunnucks, and A. C. Taylor. 2003. Demographic monitoring
of an entire species (the northern hairy-nosed wombat, Lasiorhinus kretii) by genetic analysis of
non-invasively collected material. Animal Conservation 6:101-107. DOI: 10.1017/S1367943003003135
}
\examples{
data(wombat);
wombat
}
