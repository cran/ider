#' Higher-order Local Information Dimension Estimator.
#' 
#' \code{side} is a Higher-order Information Dimension Estimator, which estimates
#' local information dimension of given dataset based on the polynomial regression
#' with Poisson error structure.
#'   
#' 
#' A variant of fractal dimension called the local information dimension is considered.
#' The local information dimension is estimated by using the probability mass function.
#' The function \code{side} considers higher-order expansion of the probability mass around
#' the inspection point, and it estimates the local information dimension by fitting 
#' a generalized linear model with Poisson error structure and an identity link function.
#' There are two methods for dimension estimation: the first method tries different 
#' dimensions and adopt the one with maximum likelihood, while the second method directly
#' maximises the likelihood with respect to the intrinsic dimension. The former returns 
#' an integer-valued dimension estimate, and the latter returns a real-valued estimate.
#' The result of the former method is used as an initial value for the latter method
#' in numerical optimization. Slow but more accurate than \code{mada} in some cases.
#'
#'
#' @inheritParams lbmle
#' @inheritParams mada
#' @param method algorithm to estimate intrinsic dimension. 'disc' for discrite dimension estimation. 'cont' for continuous dimension estimation with MLE by Newton-method.
#' @return Estimated local or global intrinsic dimension.
#' @author Hideitsu Hino \email{hideitsu.hino@@gmail.com}
#' @references H. Hino, J. Fujiki, S. Akaho, and N. Murata, 'Local Intrinsic Dimension Estimation by Generalized Linear Modeling', Neural Computation, 2017
#' @examples 
#' ## local intrinsic dimension estimate
#' tmp <- gendata(DataName='ldbl', n=300)
#' x <- tmp$x
#' set.seed(999)
#' idx <- c(sample(which(tmp$tDim==1)[1:10],3), sample(which(tmp$tDim==2)[1:30],3))
#' estmada <- mada(x=x[1:100,], local=TRUE)
#' estmada[idx]  ## estimated local intrinsic dimensions by mada
#' tmp$tDim[idx] ## true local intrinsic dimensions
#' estside <- side(x=x[1:100,], local=TRUE)
#' estside[idx] ## estimated local intrinsic dimensions by side
#' @export
side <- function(x, maxDim = 5, DM = FALSE, local = FALSE, method = "disc", comb = "average") {
    if (is.null(x)) {
        stop("data or distance matrix x is missing")
    }
    
    if (DM == FALSE) {
        x <- x/max(abs(x))
        Dmat <- as.matrix(stats::dist(x))
    } else {
        Dmat <- as.matrix(x)
    }
    
    N <- dim(Dmat)[1]
    rm(x)
    
    m <- min(round(N/5), 100)  ## number of samples used for fitting
    
    
    ## use glm2, which offers better convergence property for GLM with a non-canonical link
    ## function
    glm <- function(..., method = glm2::glm.fit2) {
        stats::glm(..., method = method)
    }
    
    
    
    estDims <- estDimsLk <- estDimsCont <- numeric(N)
    for (i in 1:N) {
        dst <- Dmat[i, ]
        dst <- dst[-i]
        Model <- list()
        ## ID <- order(dst)[1:m] ## use closest m points eps <- as.numeric(dst[ID])
        eps <- as.numeric(sort(dst, partial = m)[1:m])
        tmp <- apply(as.matrix(dst), 1, FUN = function(x) {
            x < eps
        })
        n.eps <- rowSums(tmp)
        
        if (sum(is.na(n.eps))) {
            rm.id <- which(is.na(n.eps))
            eps <- eps[-rm.id]
            tmp <- tmp[-rm.id, ]
            n.eps <- rowSums(tmp)
        }
        
        if (sum(n.eps == 0) != 0) {
            rm.id <- which(n.eps == 0)
            eps <- eps[-rm.id]
            tmp <- tmp[-rm.id, ]
            n.eps <- rowSums(tmp)
        }
        my.m <- length(eps)
        ys <- n.eps
        
        normeps <- eps/max(eps + 0.001)
        d <- 0
        
        for (my.dim in 1:maxDim) {
            d <- d + 1
            x1 <- normeps^my.dim
            x2 <- normeps^(my.dim + 2)
            regDat <- data.frame(y = ys, x1 = x1, x2 = x2)
            
            ## MLE
            model <- try(suppressWarnings(stats::glm(y ~ x1 + x2 - 1, data = regDat, fam = stats::poisson(link = "identity"), 
                start = c(0.1, 0.1))))
            
            if (class(model)[1] == "try-error") {
                stop("glm error")
            }
            Model[[length(Model) + 1]] <- model
        }
        
        estDimsLk[i] <- mID <- which.max(unlist(lapply(Model, FUN = function(x) {
            as.numeric(stats::logLik(x))
        })))
        init <- c(estDimsLk[i], stats::coef(Model[[mID]]))
        
        if (method == "cont") {
            ## (minus of) likelihood function
            lkcont <- function(arg) {
                estdim <- arg[1]
                beta0 <- arg[2]
                beta1 <- arg[3]
                x1 <- eps^estdim
                x2 <- eps^(estdim + 2)
                x <- x1 * beta0 + x2 * beta1
                -sum(-x + ys * log(x) - log(factorial(ys)))
            }
            ## (minus of) derivative of likelihood function
            dlkcont <- function(arg) {
                estdim <- arg[1]
                beta0 <- arg[2]
                beta1 <- arg[3]
                x1 <- eps^estdim
                x2 <- eps^(estdim + 2)
                x <- x1 * beta0 + x2 * beta1
                -c(sum((-x + ys) * log(eps)), sum(-x1 + ys/(beta0 + beta1 * eps^2)), sum((eps^2) * 
                  (-x1 + ys/(beta0 + beta1 * eps^2))))
            }
            res <- try(suppressWarnings(stats::optim(init, lkcont, dlkcont, method = "L-BFGS-B", 
                lower = c(1, 0, -Inf), upper = c(maxDim, Inf, Inf))), silent = TRUE)
            if (class(res) == "try-error") {
                res <- try(suppressWarnings(stats::optim(init, lkcont, dlkcont, method = "BFGS")), 
                  silent = TRUE)
                if (class(res) == "try-error") {
                  estDimsCont[i] <- init[1]
                } else if (res$par[1] > 0) {
                  estDimsCont[i] <- res$par[1]
                } else {
                  estDimsCont[i] <- init[1]
                }
            } else {
                estDimsCont[i] <- res$par[1]
            }
            
            
            estDims[i] <- estDimsCont[i]
        } else {
            estDims[i] <- estDimsLk[i]
        }
        
    }
    
    if (local == TRUE) {
        return(estDims)
    } else {
        if (comb == "average") {
            est <- mean(estDims)
        } else if (comb == "median") {
            est <- stats::median(estDims)
        } else if (comb == "vote") {
            est <- as.numeric(names(which.max(table(estDims))))
        }
        return(est)
    }
}
