#' Intrinsic Dimensionality Estimation from Near-Neighbor Information.
#' 
#' \code{nni} estimates intrinsic dimension of given dataset based on the nearest-neighbor
#' information.
#' 
#' First order expansion of the probability mass function is considered, then 
#' the distribution of the nearest-neighbor points from the inspection point is 
#' modeled by the Poisson distribution. The average of the nearest-distance is 
#' expressed by intrinsic dimension to be estimated.
#'
#'
#' @inheritParams lbmle
#' @param eps accuracy parameter. 
#' @return Estimated global intrinsic dimension.
#' @author Hideitsu Hino \email{hideitsu.hino@@gmail.com}
#' @references B. Kegl. Intrinsic dimension estimation using packing numbers.
#'  Advances in Neural Information Processing Systems 15, 2002.
#' @references K. W. Pettis et al. An intrinsic dimensionality estimator from near
#' neighbor information. IEEE transactions on pattern recognition and machine intelligence, 1979.
#' @examples
#' x <- gendata(DataName='SwissRoll',n=300)
#' estnni <- nni(x=x)
#' print(estnni)
#' @export
nni <- function(x, k1 = 2, k2 = 30, DM = FALSE, eps = 0.01, p = NULL) {
    if (is.null(x)) {
        stop("data or distance matrix x is missing")
    }
    if (DM == FALSE) {
        TkX <- function(x, k) {
            FNN::get.knn(x, k = k)$nn.dist
        }
        TkXval <- TkX(x, k2)[, k1:k2]
    } else {
        x <- as.matrix(x)
        n <- dim(x)[1]
        knn.mat <- matrix(0, ncol = k2, nrow = n)
        TkXval <- knn.mat
        for (i in 1:n) {
            knn.mat[i, ] <- order(x[i, ])[2:(k2 + 1)]
            TkXval[i, ] <- x[i, knn.mat[i, ]]
        }
        TkXval <- TkXval[, k1:k2]
    }
    RK <- colMeans(TkXval)
    
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
    
    klist <- k1:k2
    
    regDat <- data.frame(x = log(klist), y = log(RK))
    
    qestlist <- qest <- as.numeric(1/(stats::lm(y ~ 1 + x, data = regDat)$coef[2]))
    for (k in 1:100) {
        lGkd <- (qest - 1)/(2 * klist * qest^2) + (qest - 1) * (qest - 2)/(12 * (klist^2) * 
            qest^3) - ((qest - 1)^2)/(12 * (klist^3) * qest^4) - (qest - 1) * (qest - 2) * 
            (qest^2 + 3 * qest - 3)/(120 * (klist^4) * (qest^5))
        
        regDat <- data.frame(x = log(klist), y = log(RK) + lGkd)
        qest <- as.numeric(1/(stats::lm(y ~ 1 + x, data = regDat)$coef[2]))
        if (abs(qestlist[length(qestlist)] - qest) < eps) {
            return(qest)
        }
        qestlist <- c(qestlist, qest)
    }
    return(1)
    print("ERROR:did not converged")
}
