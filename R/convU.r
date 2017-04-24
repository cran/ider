#' Intrinsic Dimension Estimation with Convergence Property of a U-statistics.
#' 
#' \code{convU} estimates intrinsic dimension of given dataset based on 
#' the convergence property of Ustatistics(smoothed correlation dimension)
#'  w.r.t. kernel bandwidth
#'   
#' 
#' A variant of fractal dimension called the correlation dimension is considered.
#' The correlation dimension is defined by the notion of the correlation integral, which 
#' is calculated by counting the number of pairs closer than certain threshold epsilon.
#' The counting operation is replaced with the kernel smoothed version, and based on
#' the convergence property of the resulting U-statistics, an intrinsic dimension estimator is derived.
#'
#'
#' @inheritParams lbmle
#' @param maxDim maximum of the candidate dimension.
#' @return Estimated global intrinsic dimension.
#' @author Hideitsu Hino \email{hideitsu.hino@@gmail.com}
#' @references M. Hein and J-Y. Audibert. Intrinsic dimensionality estimation of
#' submanifolds in Rd. International Conference on Machine Learning, 2005.
#' @examples 
#' x <- gendata(DataName='SwissRoll',n=300)
#' estconvU <- convU(x=x)
#' print(estconvU)
#' @export
convU <- function(x, maxDim = 5, DM = FALSE) {
    if (is.null(x)) {
        stop("data or distance matrix x is missing")
    }
    
    if (DM == FALSE) {
        distmat <- as.matrix(stats::dist(x))
    } else {
        distmat <- as.matrix(x)
    }
    
    N <- dim(distmat)[1]
    
    oneNN <- as.numeric(apply(distmat, 2, sort)[2, ])
    
    hlN <- mean(oneNN)
    
    kh <- function(x, h) {
        ifelse((1 - (x/h)^2) > 0, 1 - (x/h)^2, 0)/h^l
    }
    slopes <- NULL
    for (l in 1:maxDim) {
        ## n=N
        Uhk1 <- mean(kh(as.numeric(stats::as.dist(distmat)), 1))
        ## n=[N/2]
        n <- round(N/2)
        hln <- hlN * ((N/n) * (log(n)/log(N)))^l
        id1 <- sample(1:N, n, replace = FALSE)
        id2 <- setdiff(1:N, id1)
        tmpD <- distmat[id1, id1]
        Uhk2 <- mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln))
        tmpD <- distmat[id2, id2]
        Uhk2 <- c(Uhk2, mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln)))
        tmpD <- distmat[id1, id2]
        Uhk2 <- c(Uhk2, mean(kh(as.numeric((tmpD)), h = hln)))
        
        Uhk2 <- mean(Uhk2)
        
        ## n=[N/3]
        n <- round(N/3)
        hln <- hlN * ((N/n) * (log(n)/log(N)))^l
        id1 <- sample(1:N, n, replace = FALSE)
        id2 <- sample(setdiff(1:N, id1), n, replace = FALSE)
        id3 <- setdiff(1:N, c(id1, id2))
        tmpD <- distmat[id1, id1]
        Uhk3 <- mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln))
        tmpD <- distmat[id2, id2]
        Uhk3 <- c(Uhk3, mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln)))
        tmpD <- distmat[id3, id3]
        Uhk3 <- c(Uhk3, mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln)))
        tmpD <- distmat[id1, id2]
        Uhk3 <- c(Uhk3, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id1, id3]
        Uhk3 <- c(Uhk3, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id2, id3]
        Uhk3 <- c(Uhk3, mean(kh(as.numeric((tmpD)), h = hln)))
        
        Uhk3 <- mean(Uhk3)
        
        ## n=[N/4]
        n <- round(N/4)
        hln <- hlN * ((N/n) * (log(n)/log(N)))^l
        id1 <- sample(1:N, n, replace = FALSE)
        id2 <- sample(setdiff(1:N, id1), n, replace = FALSE)
        id3 <- sample(setdiff(1:N, c(id1, id2)), n, replace = FALSE)
        id4 <- setdiff(1:N, c(id1, id2, id3))
        tmpD <- distmat[id1, id1]
        Uhk4 <- mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln))
        tmpD <- distmat[id2, id2]
        Uhk4 <- c(Uhk4, mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln)))
        tmpD <- distmat[id3, id3]
        Uhk4 <- c(Uhk4, mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln)))
        tmpD <- distmat[id4, id4]
        Uhk4 <- c(Uhk4, mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln)))
        
        tmpD <- distmat[id1, id2]
        Uhk4 <- c(Uhk4, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id1, id3]
        Uhk4 <- c(Uhk4, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id1, id4]
        Uhk4 <- c(Uhk4, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id2, id3]
        Uhk4 <- c(Uhk4, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id2, id4]
        Uhk4 <- c(Uhk4, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id3, id4]
        Uhk4 <- c(Uhk4, mean(kh(as.numeric((tmpD)), h = hln)))
        
        Uhk4 <- mean(Uhk4)
        
        ## n=[N/5]
        n <- round(N/5)
        hln <- hlN * ((N/n) * (log(n)/log(N)))^l
        id1 <- sample(1:N, n, replace = FALSE)
        id2 <- sample(setdiff(1:N, id1), n, replace = FALSE)
        id3 <- sample(setdiff(1:N, c(id1, id2)), n, replace = FALSE)
        id4 <- sample(setdiff(1:N, c(id1, id2, id3)), n, replace = FALSE)
        id5 <- setdiff(1:N, c(id1, id2, id3, id4))
        tmpD <- distmat[id1, id1]
        Uhk5 <- mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln))
        tmpD <- distmat[id2, id2]
        Uhk5 <- mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln))
        tmpD <- distmat[id3, id3]
        Uhk5 <- mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln))
        tmpD <- distmat[id4, id4]
        Uhk5 <- mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln))
        tmpD <- distmat[id5, id5]
        Uhk5 <- mean(kh(as.numeric(stats::as.dist(tmpD)), h = hln))
        
        tmpD <- distmat[id1, id2]
        Uhk5 <- c(Uhk5, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id1, id3]
        Uhk5 <- c(Uhk5, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id1, id4]
        Uhk5 <- c(Uhk5, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id1, id5]
        Uhk5 <- c(Uhk5, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id2, id3]
        Uhk5 <- c(Uhk5, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id2, id4]
        Uhk5 <- c(Uhk5, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id2, id5]
        Uhk5 <- c(Uhk5, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id3, id4]
        Uhk5 <- c(Uhk5, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id3, id5]
        Uhk5 <- c(Uhk5, mean(kh(as.numeric((tmpD)), h = hln)))
        tmpD <- distmat[id4, id5]
        Uhk5 <- c(Uhk5, mean(kh(as.numeric((tmpD)), h = hln)))
        
        Uhk5 <- mean(Uhk5)
        
        ns <- round(N/1:5)
        x <- log(hlN * ((N/ns) * (log(ns)/log(N)))^l)
        y <- log(c(Uhk1, Uhk2, Uhk3, Uhk4, Uhk5))
        
        xid <- is.finite(x)
        yid <- is.finite(y)
        x <- x[which((xid * yid) == 1)]
        y <- y[which((xid * yid) == 1)]
        
        regDat <- data.frame(x = x, y = y)
        tmp <- try(stats::lm(y ~ x + 1, data = regDat, weight = (1/(1:5))[which((xid * yid) == 
            1)])$coef[2])
        slopes <- c(slopes, tmp)
    }
    
    return(max(which.min(abs(slopes)), 1))
}
