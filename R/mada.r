#'  Manifold-Adaptive Local Dimension Estimation.
#' 
#' \code{mada} estimates local information dimension of given dataset based on 
#' the first order expansion of probability mass function.
#'  
#' 
#' A variant of fractal dimension called the local information dimension is considered.
#' The local information dimension is estimated by using the probability mass function.
#' The function \code{mada} considers first order expansion of the probability mass around
#' the inspection point, and it estimates the local information dimension by using two different
#' radii from the inspection point.
#'
#'
#' @inheritParams lbmle
#' @param k k-NN parameter.
#' @param comb 'average', 'median' or 'vote' for combining local estimates when global estimate is required.
#' @param local logical. If \code{TRUE}, a vector of local dimensions at each sample point is returned.
#' @param maxDim maximum of the candidate dimensions.
#' @return Estimated local or global intrinsic dimension.
#' @author Hideitsu Hino \email{hideitsu.hino@@gmail.com}
#' @references A. M. Farahmand, C. Szepesvari and J-Y. Audibert.
#'  Manifold-adaptive dimension estimation. International Conference on Machine Learning, 2007.
#' @examples 
#' ## local intrinsic dimension estimate
#' tmp <- gendata(DataName='ldbl',n=300)
#' x <- tmp$x
#' estmada <- mada(x=x,local=TRUE)
#' head(estmada)  ## estimated local intrinsic dimensions by mada
#' head(tmp$tDim) ## true local intrinsic dimensions
#' @export
mada <- function(x, k = NULL, comb = "average", DM = FALSE, local = FALSE, maxDim = 5) {
    if (is.null(x)) {
        stop("data or distance matrix x is missing")
    }
    
    
    if (DM == FALSE) {
        distmat <- as.matrix(stats::dist(x))
    } else {
        distmat <- as.matrix(x)
    }
    n <- dim(distmat)[1]
    
    if (is.null(k)) {
        k <- floor(2 * log(n))
    }
    
    if (local == FALSE) {
        ## sample m=round(n/2) point as centers, and estimate the local dimensionality to be
        ## combined
        id <- sample(1:n, round(n/2), replace = FALSE)
        tmpD <- distmat[id, ]
        tmpD[which(tmpD == 0)] <- max(tmpD)
    } else {
        tmpD <- distmat
        tmpD[which(tmpD == 0)] <- max(tmpD)
    }
    
    ## for sampled points, compute k-th and [k/2]-th nearest distances
    sortedD <- apply(tmpD, 1, sort)
    RK <- sortedD[k, ]
    RK2 <- sortedD[floor(k/2), ]
    
    ests <- log(2)/log(RK/RK2)
    
    ## return estimated dimensions of each point
    if (local == TRUE) {
        ## return(round(ests))
        return(ests)
    }
    
    if (comb == "average") {
        est <- mean(ests)
    } else if (comb == "median") {
        est <- stats::median(ests)
    } else {
        est <- as.numeric(names(which.max(table(ests))))
    }
    return(est)
}
