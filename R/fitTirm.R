#' @title Fit Two Innate Rates Model (TIRM)
#'
#' @description Fits the Two Innate Rates Model (TIRM) to count data to obtain the MLE for population size
#'
#' @param data A two-column data frame with the first column specifiying the capture class (i.e. individuals in class i were caught i times) and the second column specifying the number of individuals in each class
#'
#' @param max.pop The maximum population size
#'
#' @param max.iter The maximum number of iterations to acheive convergence. The default will be more than enough for most cases
#'
#' @details The TIRM model fit by this function assumes that individuals can be assigned to two classes. Class A represent the frequently captured individuals. Class B represents the infrequently captured individuals.
#'
#' The value is specified for \code{max.pop} is not likely to matter as long as it is much greater than the maximum likelihood estimate for population size.
#'
#' Note that if the data contains only singletons, the data is not informative and the maximum likelihood estimate for population size will be equal to \code{max.pop}
#'
#' @return an object of class \code{capfit} containing the follow elements:
#' \describe{
#'  \item{model}{The model specified}
#'  \item{likelihood}{The likelihood of the model}
#'  \item{ml.pop.size}{The maximum likelihood estimate for population size}
#'  \item{ml.na}{The maximum likelihood estimate for the number of individuals in class A}
#'  \item{ml.nb}{The maximum likelihood estimate for the number of individuals in class B}
#'  \item{alpha}{The ratio of the rates of captures between class A and class B individuals}
#'  \item{cap.ind}{The mean number of captures per individual}
#'  \item{sampled.ind}{The total number of individuals in the sample}
#'  \item{sample.size}{Total number of samples in the data set}
#'  \item{max.pop}{The maximum population size specified by \code{max.pop}}
#' }
#'
#' @references Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005.
#'
#' Pennell, M. W., C. R. Stansbury, L. P. Waits and C. R. Miller. 2013. Capwire: a R package for estimating population census size from non-invasive genetic sampling. Molecular Ecology Resources 13:154-157. 
#'
#' @seealso \code{\link{simTirm}}
#'
#' @author Matthew W. Pennell
#'
#' @export fitTirm
#'
#' @examples
#' ## Simulate data under Two Innate Rates Model
#'
#' data <- simTirm(na=20, nb=15, a=4, s=150)
#'
#' ## Fit Two Innate Rates Model to Data
#'
#' res <- fitTirm(data=data, max.pop=200)
#' res
#' 
fitTirm <- function(data, max.pop, max.iter=20){
    check.capwire.data(data)

    ## get basic info from data
    si <- get.sampling.info(data)

    ## return object if only singletons
    if (all(si$counts == 1)){
        r <- list(model="Two.innate.rates", likelihood=NA,
                  ml.pop.size=max.pop, ml.na=0, ml.nb=si$t,
                  alpha=NA, cap.ind=si$s/si$t, sampled.ind=si$t,
                  sample.size=si$s, max.pop=max.pop)
        class(r) <- c("list", "capfit")
        return(r)
    }
    
    ## initialize alpha (step 1) 
    a <- init.alpha(si)
   
    ## set initial variables
    na <- round(si$t/2)
    nb <- si$t - na
    cc1 <- list(ca=0, cb=c(0,0))
    cc2 <- list(ca=1, cb=c(1))
    i <- 0

    ## iterate
    while(1){
        if (length(cc1$cb) == length(cc2$cb) | i == max.iter)
            break()

        ## assign capture classes (step 3)
        cc1 <- capture.classes(na, nb, si, a)
        ## calculate likelihood (step 4)
        ml <- find.mle.tirm(cc1, a, max.pop)
        ## correct bias (step 5)
        a <- adjust.alpha(ml, si$s)
        ## calculate likelihood (step 6)
        ml <- find.mle.tirm(cc1, a, max.pop)
        ## reassign (step 7) 
        cc2 <- capture.classes(ml$na, ml$nb, si, a)

        ## reset variables
        na <- ml$na
        nb <- ml$nb
        i <- i + 1
        
    }

    ## add constant to likelihood
    ## need to check this
    lik <- ml$lik - sum(log(seq_len(length(cc2$ca)))) -
        sum(log(seq_len(length(cc2$cb))))
    
    r <- list(model="Two.innate.rates", likelihood=lik,
              ml.pop.size=(ml$na+ml$nb), ml.na=ml$na, ml.nb=ml$nb,
              alpha=a, cap.ind=si$s/si$t, sampled.ind=si$t,
              sample.size=si$s, max.pop=max.pop)
    class(r) <- c("list", "capfit")
    r
}
        



## initialize alpha
## takes object of class sampling info
init.alpha <- function(si){
    cap.ind <- si$s/si$t
    ab <- mean(si$counts[si$counts > cap.ind])
    bl <- mean(si$counts[si$counts <= cap.ind])
    ab/bl
}


## get expected capture counts
expected.capture.counts <- function(na, nb, s, a){
    eca <- s * a /(a * na + nb)
    ecb <- s / (a * na + nb)
    list(eca=eca, ecb=ecb)
}

## assign to capture classes
capture.classes <- function(na, nb, si, a){
    counts <- si$counts
    cb <- counts[counts == 1]
    rem <- counts[counts != 1]
    ec <- expected.capture.counts(na, nb, si$s, a)
    ## assign all other individuals to one of the classes
    asg <- sapply(rem, function(x) {abs(x - ec$eca) < abs(x - ec$ecb)})
    ca <- rem[asg]
    cb <- c(cb, rem[!asg])
    list(ca=ca, cb=cb)
}


## adjust the bias in alpha
adjust.alpha <- function(fit, s)
    (fit$nb * sum(fit$ca)) / (fit$na * s - fit$na * sum(fit$ca))


## find mle
find.mle.tirm <- function(cc, a, max.pop){
    ta <- na <- length(cc$ca)
    tb <- nb <- length(cc$cb)
    v <- c(nb:(max.pop - na))
    lik <- sapply(v, function(x)
                  lik.tirm(na, x, a, cc$ca, cc$cb, ta, tb))
    nb <- v[which(lik == max(lik))]
    list(lik=max(lik), na=na, nb=nb, ca=cc$ca, cb=cc$cb)
}


## likelihood function for tirm
lik.tirm <- function(na, nb, a, ca, cb, ta, tb){
    lik <- sum(log(seq_len(na))) + log(a/(a*na + nb)) * sum(ca) +
            sum(log(seq_len(nb))) + log(1/(a*na + nb)) *sum(cb)

    if (nb != tb)
        lik <- lik - sum(log(seq_len(nb - tb)))

    lik
}

    
            
    
    


