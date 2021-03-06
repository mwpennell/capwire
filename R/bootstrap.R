#' @title Parametric bootstrap ECM and TIRM models
#'
#' @description Uses maximum likelihood parameter estimates from \code{fitEcm}, or \code{fitTirm} to perform a parametric bootstrap to get confidence intervals for the estimate of the population size
#'
#' @param fit a \code{capfit} object inherited from \code{\link{fitEcm}} or \code{\link{fitTirm}}
#'
#' @param bootstraps The number of bootstraps to be performed (default is 1000)
#'
#' @param CI A vector of quantiles to to generate a confidence interval for the population estimate. The default is \code{c(0.025, 0.975)}, denoting a symmetrical 95 percent confidence interval.
#'
#' @details This function uses the ML estimates obtained from fitting the model to simulate data under the model.
#'
#' \code{\link{bootstrapCapwire}} inherits an object from \code{\link{fitEcm}} or \code{\link{fitTirm}} such that the model and parameter estimates do not need to be specified.
#' 
#' The ML estimate for the population size will also be returned but this will not be changed by \code{bootstrapCapwire}
#' 
#' Note that if the model is a poor fit to the data, the confidence intervals may not be reliable.
#' 
#' The lower confidence interval is bounded by the number of unique individuals in the sample.
#'
#' @return A list containing the following elements:
#' \describe{
#'  \item{ml.pop.size}{The maximum likelihood estimate for the population size obtained by fitting the model}
#'  \item{conf.int}{The confidence interval for the estimate of the population size}
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
#' @export bootstrapCapwire
#'
#' @examples
#' ## Simulate data under Equal Capture Model
#'
#' data <- simEcm(n=40, s=150)
#'
#' ## Fit Equal Capture Model to Data
#'
#' res <- fitEcm(data=data, max.pop=200)
#'
#' ## Perform bootstrap to get confidence intervals
#'
#' ci <- bootstrapCapwire(fit=res, bootstraps=50, CI = c(0.025, 0.975))
#' ci
#' 
bootstrapCapwire <- function(fit, bootstraps=1000, CI=c(0.025, 0.975)){
    if (!inherits(fit, "capfit"))
        stop("'fit' needs to be a capwire model fitted object")

    model <- fit$model
    boot <- switch(model,
                   Equal.capture=bootstrap.ecm(fit, bootstraps),
                   Two.innate.rates=bootstrap.tirm(fit, bootstraps))

    ci <- quantile(boot, CI)
    list(ml.pop.size=fit$ml.pop.size, conf.int=ci)
}



bootstrap.ecm <- function(fit, bootstraps){
    n <- fit$ml.pop.size
    s <- fit$sample.size
    sims <- lapply(seq_len(bootstraps), function(x) simEcm(n, s))
    fits <- lapply(sims, function(x) {fitEcm(x, fit$max.pop)$ml.pop.size})
    do.call(c,fits)
}

bootstrap.tirm <- function(fit, bootstraps){
    na <- fit$ml.na
    nb <- fit$ml.nb
    a <- fit$alpha
    s <- fit$sample.size
    sims <- lapply(seq_len(bootstraps), function(x) simTirm(na, nb, a, s))
    fits <- lapply(sims, function(x) {fitTirm(x, fit$max.pop)$ml.pop.size})
    do.call(c,fits)
}



