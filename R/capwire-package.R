#' Capwire: estimate population census size from non-invasive genetic sampling
#'
#' @description This R package is designed to fit the sampling models developed by Miller et al. 2005 (Molecular Ecology). The two main model fitting functions are \code{\link{fitTirm}} and \code{\link{fitEcm}}. The package also includes tools to estimate confidence intervals by parameteric bootstrapping (\code{\link{bootstrapCapwire}}), perform model selection (\code{\link{lrtCapwire}}) and simulate data under a variety of conditions (\code{\link{simTirm}}, \code{\link{simEcm}}, and \code{\link{simCapture}}).
#'
#' @references Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005.
#'
#' Pennell, M. W., C. R. Stansbury, L. P. Waits and C. R. Miller. 2013. Capwire: a R package for estimating population census size from non-invasive genetic sampling. Molecular Ecology Resources 13:154-157.
#'
#' @docType package
#' @name capwire
#' @aliases capwire-capwire capwire
#'
NULL

#' @title Wombat dataset
#'
#' This dataset is used to demonstrate the functionality of the package (see supplemental material of Pennell et al. 2013). It can be used as input for \code{\link{fitTirm}} and \code{\link{fitEcm}}.
#'
#' @references Banks, S. C., S. D. Hoyle, A. Horsup, P. Sunnucks, and A. C. Taylor. 2003. Demographic monitoring of an entire species (the northern hairy-nosed wombat, Lasiorhinus kretii) by genetic analysis of non-invasively collected material. Animal Conservation 6:101-107.
#'
#' @docType data
#' @keywords datasets
#' @format data frame
#' @name wombat
NULL
