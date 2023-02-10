#' Maximum Likelihood Estimation of Intrinsic Dimension.
#'
#' \code{lbmle} estimate the intrinsic dimension of a given dataset.
#' 
#' The likelihood of the rate parameter of the Poisson process, which characterize the behaviour of  
#' the distance from a point to another point in the given dataset, is considered, and the maximum likelihood estimator (MLE) for the intrinsic dimension is derived.
#' The original method proposed by Levina and Bickel contains a known bias, and it is corrected by Mackay and Ghahramani. This function implements both, with the default the bias corrected estimate.
#'
#'
#' @param x data matrix or distance matrix given by as.matrix(dist(x)).
#' @param k1 first k-NN parameter. 
#' @param k2 second k-NN parameter.
#' @param BC whether bias is corrected or not. logical.
#' @param DM whether \code{'x'} is distance matrix or not. logical.
#' @param p ambient dimension used for automatically define \code{'k1'} and \code{'k2'}.
#' @return Estimated global intrinsic dimension.
#' @author Hideitsu Hino \email{hideitsu.hino@@gmail.com}
#' @references E. Levina and P. J. Bickel. Maximum likelihood estimation of 
#' intrinsic dimension. Advances in Neural Information Processing Systems 17, 2005.
#' @references D. MacKay and Z. Ghahramani. \url{http://www.inference.org.uk/mackay/dimension/}
#' @examples 
#' x <- gendata(DataName='SwissRoll',n=300)
#' estmle <- lbmle(x=x,k1=3,k2=5)
#' print(estmle)
#' @export
lbmle <- function(x = NULL, k1 = NULL, k2 = NULL, BC = TRUE, DM = FALSE, p = NULL) {
    if (is.null(x)) {
        stop("data or distance matrix x is missing")
    }
    if (DM == TRUE) {
        if (inherits(x,"dist")) {
            x <- as.matrix(x)
        }
    }
    n <- dim(x)[1]  ## number of observations
    
    ## if k is not given, we use a rule of thumb from the asymptotic optimal order
    if (is.null(k1) && !is.null(p)) {
        k1 <- round(0.5 * n^(4/(4 + p)))
    }
    if (is.null(k1) && is.null(p)) {
        stop("specify the ambient dimension p, if k1 and k2 are omitted")
    }
    if (is.null(k2) && !is.null(p)) {
        k2 <- round(n^(4/(4 + p)))
    }
    if (is.null(k2) && is.null(p)) {
        stop("specify the ambient dimension p")
    }
    
    if (DM == FALSE) {
        TkX <- function(x, k) {
            FNN::get.knn(x, k = k)$nn.dist
        }
        
        hat_qk <- NULL
        for (ki in k1:k2) {
            TkXval <- TkX(x, ki)
            if (BC) {
                hat_qk <- c(hat_qk, 1/mean(log(TkXval[, ki]) - (rowSums(log(TkXval)[, 1:(ki - 
                  1)])/(ki - 1))))
            } else {
                hat_qk <- c(hat_qk, mean(1/(log(TkXval[, ki]) - rowSums(log(TkXval)[, 1:(ki - 
                  1)])/(ki - 1))))
            }
        }
    } else {
        
        hat_qk <- NULL
        for (ki in k1:k2) {
            knn.mat <- matrix(0, ncol = ki, nrow = n)
            TkXval <- knn.mat
            for (i in 1:n) {
                knn.mat[i, ] <- order(x[i, ])[2:(ki + 1)]
                TkXval[i, ] <- x[i, knn.mat[i, ]]
            }
            if (BC) {
                hat_qk <- c(hat_qk, 1/mean(log(TkXval[, ki]) - (rowSums(log(TkXval)[, 1:(ki - 
                  1)])/(ki - 1))))
            } else {
                hat_qk <- c(hat_qk, mean(1/(log(TkXval[, ki]) - rowSums(log(TkXval)[, 1:(ki - 
                  1)])/(ki - 1))))
            }
        }
        
    }
    hat_q <- sum(hat_qk)/(k2 - k1 + 1)
    return(hat_q)
}
