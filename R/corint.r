#' Intrinsic Dimension Estimation with Correlation Integral
#' 
#' \code{corint} estimates intrinsic dimension of given dataset based on 
#' the correlation integral
#'   
#' 
#' A variant of fractal dimension called the correlation dimension is considered.
#' The correlation dimension is defined by the notion of the correlation integral, 
#' is calculated by using the power low for the definition of the correlation dimension.
#'
#'
#' @inheritParams lbmle
#' @return Estimated global intrinsic dimension.
#' @author Hideitsu Hino \email{hideitsu.hino@@gmail.com}
#' @references P. Grassberger and I. Procaccia. Measuring the strangeness of strange attractors. 
#' Physica, 1983.
#' @examples 
#' x <- gendata(DataName='SwissRoll',n=300)
#' estcorint <- corint(x=x,k1=5,k2=10)
#' print(estcorint)
#' @export
corint <- function(x, k1 = NULL, k2 = NULL, DM = FALSE, p = NULL) {
    if (is.null(x)) {
        stop("data or distance matrix x is missing")
    }
    
    if (DM == FALSE) {
        distmat <- as.matrix(stats::dist(x))
    } else {
        distmat <- as.matrix(x)
    }
    n <- dim(distmat)[1]
    
    ## if k is not given, we use a rule of thumb from the asymptotic optimal order
    if (is.null(k1) && !is.null(p)) {
        k1 <- round(0.5 * n^(4/(4 + p)))
    }
    if (is.null(k1) && is.null(p)) {
        stop("specify the ambient dimension p")
    }
    if (is.null(k2) && !is.null(p)) {
        k2 <- round(n^(4/(4 + p)))
    }
    if (is.null(k2) && is.null(p)) {
        stop("specify the ambient dimension p")
    }
    
    
    ## calculate r1 and r2 from k1 and k2
    tmp <- as.vector(stats::as.dist(distmat))
    knn.mat <- matrix(0, ncol = k2, nrow = n)
    TkXval <- knn.mat
    for (i in 1:n) {
        knn.mat[i, ] <- order(distmat[i, ])[2:(k2 + 1)]
        TkXval[i, ] <- distmat[i, knn.mat[i, ]]
    }
    r1 <- stats::median(TkXval[, k1])
    r2 <- stats::median(TkXval[, k2])
    Cr <- c(mean(tmp < r1), mean(tmp < r2))
    estq <- diff(log(Cr))/log(r2/r1)
    
    return(estq)
}
