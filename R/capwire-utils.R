

#' @title Convert capture counts to table of capture classes
#'
#' @description Converts a vector of capture counts into a two-column matrix consisting of all capture classes and the individuals associated with each class.
#'
#' @param counts a vector of capture count data
#'
#' @return A two-column matrix with the first column specifiying the capture class (where all individuals in class i were caught i times) and the second column specifying the number of individuals in this capture class.
#'
#' The data can be used as the data argument for any of the model-fitting functions implemented in capwire
#'
#' @references Miller C. R., P. Joyce and L.P. Waits. 2005. A new method for estimating the size of small populations from genetic mark-recapture data. Molecular Ecology 14:1991-2005.
#'
#' Pennell, M. W., C. R. Stansbury, L. P. Waits and C. R. Miller. 2013. Capwire: a R package for estimating population census size from non-invasive genetic sampling. Molecular Ecology Resources 13:154-157. 
#'
#' @seealso \code{\link{fitTirm}}, \code{\link{fitEcm}}
#'
#' @author Matthew W. Pennell
#'
#' @export buildClassTable
#'
#' @examples
#'
#' ## create a vector of capture counts
#' 
#' counts <- c(1,1,1,1,1,2,2,3,3,4,5)
#'
#' ## build table 
#' 
#' d <- buildClassTable(counts)
#' d
#'
buildClassTable <- function(counts){
    if (!inherits(counts, "numeric"))
        stop("counts needs to be a numeric vector")
    
    uni <- unique(counts)
    r <- sort(sapply(uni, function(x) length(counts[counts == x])))
    res <- lapply(unique(r), function(x) c(x, length(r[r == x])))
    res <- do.call(rbind, res)
    colnames(res) <- c("capture.class", "n.ind")
    res
}


## check capwire object
check.capwire.data <- function(x){
    if (!"matrix" %in% class(x) & !"data.frame" %in% class(x))
        stop("data must be entered as either a 'data.frame' or 'matrix'")

    if (ncol(x) != 2)
        stop("data should include exactly two columns")

}


## get sampling info
get.sampling.info <- function(d){
    counts <- lapply(seq_len(nrow(d)), function(x) {rep(d[x,1], d[x,2])})
    counts <- do.call(c, counts)
    ## remove 0 values
    counts <- counts[counts > 0]
    s <- sum(counts)
    t <- length(counts)
    list(counts=counts, s=s, t=t)
}
