#' @title Likelihood ratio test comparing TIRM to ERM
#'
#' @description Performs Likelihood Ratio Test to compare model fits between the Equal Capture Model (ECM) and Two Innate-Rates Model (TIRM).
#' 
#' A parametric bootstrap is used to generate a null-distribution for the likelihood ratio test-statistic to assess significance
#'
#' @param ecm A \code{capfit} object inherited from \code{\link{fitEcm}}
#'
#' @param tirm A \code{capfit} object inherited from \code{\link{fitTirm}}
#'
#' @param bootstraps The number of parameteric bootstraps to perform to generate the null distribution
#'
#' @details The maximum likelihood estimate for population size under the ECM model is used to simulate data under the ECM model. The data is then fit to both the ECM and TIRM models and the likelihood ratio test is evaluated.
#'
#' The observed Likelihood Ratio statistic is then compared to the distribution of likelihood ratio statistics observed when the null model is true in order to estimate a p-value.
#'
#' In most applications of likelihood ratio tests, the more parameterized model will always have a higher likelihood. However, this is not necessarily the case for the models implemented in \code{capwire}.The likelihood involves a combinatorial term which denotes the number of ways a population could give rise to the observed vector of counts. As such, for some data sets (such as those for which the assumptions of the ECM model are valid),there will be more ways to obtain the observed data if the population is not subdivided into classes.
#'
#' @return A list consisting of the following items
#' \describe{
#'  \item{LR}{Likelihood ratio}
#'  \item{p.value}{P-value for the siginficance of the likelihood ratio test derived from a parametric bootstrap}
#' }
#'
#' @references Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005.
#'
#' Pennell, M. W., C. R. Stansbury, L. P. Waits and C. R. Miller. 2013. Capwire: a R package for estimating population census size from non-invasive genetic sampling. Molecular Ecology Resources 13:154-157. 
#'
#' @seealso \code{\link{fitTirm}}, \code{\link{fitEcm}}
#'
#' @author Matthew W. Pennell
#'
#' @export lrtCapwire
#'
#' @examples
#' ## Simulate data under Equal Capture Model
#' 
#' data <- simEcm(n=20, s=60)
#' 
#' ## Fit Equal Capture Model (ECM) to count data
#'
#' ecm <- fitEcm(data, max.pop=200)
#'
#' ## Fit Two-Innate Rates Model (TIRM) to count data
#' 
#' tirm <- fitTirm(data, max.pop=200)
#'
#' ## Perform Likelihood Ratio Test
#'
#' lrtCapwire(ecm=ecm, tirm=tirm, bootstraps=10)
#'
lrtCapwire <- function(ecm, tirm, bootstraps=100){
    if (!inherits(ecm, "capfit") | !inherits(tirm, "capfit"))
        stop("objects need to be capwire model fitted object")

    lr <- tirm$likelihood - ecm$likelihood
    sims <- lapply(seq_len(bootstraps), function(x)
                   simEcm(ecm$ml.pop.size, ecm$sample.size))
    lecm <- lapply(sims, function(x) fitEcm(x, ecm$max.pop))
    ltir <- lapply(sims, function(x) fitTirm(x, ecm$max.pop))
    rat <- lapply(seq_len(bootstraps), function(x)
                  {ltir[[x]]$likelihood - lecm[[x]]$likelihood})

    p <- (length(rat[rat >= lr]) + 1) / (bootstraps + 1)
    list(LR=lr, p.value=p)
}
