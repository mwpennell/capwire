#' @title Simulate data under Two Innate Rates Model (TIRM)
#'
#' @description Simulates capture count data under the assumptions of the Two-Innate Rates Model (TIRM) where the individuals in the population are assumed to come from two classes: easy to capture individuals and difficult to capture individuals.
#'
#' @param na Number of individuals in the "easy to capture" class
#'
#' @param nb Number of individuals in the "difficult to capture" class
#'
#' @param a The ratio of capture rates between the two class (i.e. the alpha parameter in \code{\link{fitTirm}}
#'
#' @param s Total number of samples collected
#'
#' @return A two-column matrix with the first column specifiying the capture class (i.e. individuals caught i times) and the second column specifying the number of individuals in each class.
#' 
#' The matrix can be used as the data argument for any of the model-fitting functions implemented in capwire
#'
#' @references Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005.
#'
#' Pennell, M. W., C. R. Stansbury, L. P. Waits and C. R. Miller. 2013. Capwire: a R package for estimating population census size from non-invasive genetic sampling. Molecular Ecology Resources 13:154-157. 
#'
#' @seealso \code{\link{fitTirm}}
#'
#' @author Matthew W. Pennell
#'
#' @export simTirm
#'
#' @examples 
#' ## Simualte data under Two-Innate Rates Model
#'
#' data <- simTirm(na=20, nb=10, a=4, s=100)
#'
#' ## Fit the Two-Innate Rates Model 
#' 
#' tirm <- fitTirm(data=data, max.pop=200)
#' tirm
#' 
simTirm <- function(na, nb, a, s){
    v <- c(seq_len(na+nb))
    d <- sample(v, size=s, prob=c(rep(1,na), rep(a, nb)), replace=TRUE)
    uni <- unique(d)
    r <- sort(sapply(uni, function(x) length(d[d == x])))
    res <- lapply(unique(r), function(x) c(x, length(r[r == x])))
    res <- do.call(rbind, res)
    colnames(res) <- c("capture.class", "n.ind")
    res
}

#' @title Simulate data under Equal Capture Model (TIRM)
#'
#' @description Simulates capture count data under the assumptions of the Equal Capture Model (ECM) where all individuals in the population are assumed to be equally likely to be captured.
#'
#' @param n Number of individuals in the population
#'
#' @param s Total number of samples collected
#'
#' @return A two-column matrix with the first column specifiying the capture class (i.e. individuals caught i times) and the second column specifying the number of individuals in each class.
#' 
#' The matrix can be used as the data argument for any of the model-fitting functions implemented in capwire
#'
#' @references Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005. Pennell, M. W., C. R. Stansbury, L. P. Waits and C. R. Miller. Capwrie: a R package for estimating population census size from non-invasive genetic sampling
#'
#' @seealso \code{\link{fitEcm}}
#'
#' @author Matthew W. Pennell
#'
#' @export simEcm
#'
#' @examples 
#' ## Simualte data under Equal Capture Model
#'
#' data <- simEcm(n=25, s=100)
#'
#' ## Fit the Equal Capture Model 
#' 
#' ecm <- fitEcm(data=data, max.pop=200)
#' ecm
#' 
simEcm <- function(n, s){
    d <- sample(seq_len(n), size=s, replace=TRUE)
    uni <- unique(d)
    r <- sort(sapply(uni, function(x) length(d[d == x])))
    res <- lapply(unique(r), function(x) c(x, length(r[r == x])))
    res <- do.call(rbind, res)
    colnames(res) <- c("capture.class", "n.ind")
    res
}





#' @title simCapture
#'
#' @aliases drawCapRatesUnif, drawCapRatesExp, drawCapRatesGamma, drawCapRatesGeom
#' 
#' @description Simulates capture count data where individual capture rates are assumed to be drawn from a specified distribution.
#' 
#' Data can be used as input for fitting Equal Capture Model (with \code{fitEcm}) or Two Innate Rates Model (with \code{fitTirm})
#'
#' @param n number of individuals in the population
#'
#' @param s total number of samples collected
#'
#' @param dist.func The distribution of capture rates within the population (see details)
#'
#' @param return.cap.probs Logical, signifying whether individual capture probabilities should be returned
#'
#' @details We assume that there is heterogeneity in the capturabilities of individuals within a population. That is, some individuals are more likely to be captured than others
#' 
#' We also assume that the individual capturabilities are drawn from some distribution.
#' 
#' The distribution is specified by the \code{dist.func} argument. \code{dist.func} takes a function with parameter n, where n specifies the number of samples to be drawn.\code{simCapture} can take any distribution of this form but the \code{capwire} package includes several functions which allow for users to draw capture rates from several standard distribution such as a uniform (\code{drawCapRatesUnif}), exponential (\code{drawCapRatesExp}), gamma (\code{drawCapRatesGamma}), geometric (\code{drawCapRatesGeom}), and beta (\code{drawCapRatesBeta}).
#'
#' @return If \code{return.cap.probs=FALSE}: A two-column matrix with the first column specifiying the capture class (i.e. individuals caught i times) and the second column specifying the number of individuals in each class.
#' 
#' If \code{return.cap.probs=TRUE}, an additional matrix is returned with the capture probabilites of every individual in the population
#'
#' #' @references Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005. Pennell, M. W., C. R. Stansbury, L. P. Waits and C. R. Miller. Capwrie: a R package for estimating population census size from non-invasive genetic sampling
#'
#' @seealso \code{\link{simTirm}}, \code{\link{simEcm}}
#'
#' @author Matthew W. Pennell
#'
#' @export simCapture
#'
#' @examples
#'
#' ## Specify a uniform distribution
#' ud <- drawCapRatesUnif(0,1)
#' simCapture(n=20, s=100, ud)
#'
#' ## Specify an exponential distribution
#' ed <- drawCapRatesExp(0.5)
#' simCapture(n=20, s=100, ed)
#'
#' ## Specify a gamma distribution
#' gd <- drawCapRatesGamma(1,0.5)
#' simCapture(n=20, s=100, gd)
#'
#' ## Specify a geometric distribution
#' md <- drawCapRatesGeom(0.5)
#' simCapture(n=20, s=100, md)
#'
#' ## Specify a beta distribution
#' bd <- drawCapRatesBeta(1, 0.5)
#' simCapture(n=20, s=100, bd)
#'
#' ## Specify a custom distribution
#' ## Here a one-tailed normal distribution
#' drawCapRatesAbsNorm <- function(mean, sd){function(n){abs(rnorm(n, mean, sd))}}
#' nd <- drawCapRatesAbsNorm(0,1)
#' simCapture(n=20, s=100, nd)
#' 
simCapture <- function(n, s, dist.func, return.cap.probs=FALSE){
    prob <- dist.func(n)
    d <- sample(seq_len(n), size=s, prob=prob, replace=TRUE)
    uni <- unique(d)
    r <- sort(sapply(uni, function(x) length(d[d == x])))
    res <- lapply(unique(r), function(x) c(x, length(r[r == x])))
    res <- do.call(rbind, res)
    colnames(res) <- c("capture.class", "n.ind")
    res

    if (return.cap.probs==TRUE)
        res <- list(data=res, cap.probs=cbind(ind=seq_len(n), cap.prob=prob))

    res
}


#' @export drawCapRatesUnif
#' @aliases simCapture
drawCapRatesUnif <-function(lower, upper){
	function(n) runif(n, min=lower, max=upper)
}

#' @export drawCapRatesGeom
#' @aliases simCapture
drawCapRatesGeom <- function(p){
	function(n) rgeom(n, prob=p)
}

#' @export drawCapRatesGamma
#' @aliases simCapture
drawCapRatesGamma <- function(shape, rate){
	function(n) rgamma(n, shape=shape, rate=rate)
}

#' @export drawCapRatesExp
#' @aliases simCapture
drawCapRatesExp <- function(r){
	function(n) rexp(n, rate=r)
}

#' @export drawCapRatesBeta
#' @aliases simCapture
drawCapRatesBeta <- function(shape1, shape2){
	function(n) rbeta(n, shape1, shape2)
}





        


    



















