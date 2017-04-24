#' Intrinsic Dimension Estimation Using Packing Numbers.
#' 
#' \code{pack} estimates intrinsic dimension of given dataset based on the packing number.
#' 
#' A variant of fractal dimension called the capacity dimension is considered.
#' The capacity dimension is defined by using the notion of covering number, which is
#' hard to calculate in general. In this function, the packing number of the data
#'  space is used as the surrogate of the covering number.
#'   The packing number is estimated by greedy manner or by hierarchical clustering.
#'  
#'
#'
#' @inheritParams lbmle
#' @param greedy logical. If \code{TRUE}, then a greedy algorithm is used for 
#' estimating the packing number. If \code{FALSE}, then a hierarchical clustering
#' algorithm is used instead.
#' @param eps accuracy parameter. 
#' @param k1 first radius parameter. If one of \code{k1} or \code{k2} is \code{NULL},
#' then both are automatically determined from the input data.
#' @param k2 second radius parameter.
#' @return Estimated global intrinsic dimension.
#' @author Hideitsu Hino \email{hideitsu.hino@@gmail.com}
#' @references B. Kegl. Intrinsic dimension estimation using packing numbers.
#'  Advances in Neural Information Processing Systems 15, 2002.
#' @references B. Eriksson and M. Crovella. Estimating intrinsic dimension via clustering.
#' IEEE Statistical Signal Processing Workshop, 2012.
#' @examples
#' x <- gendata(DataName='SwissRoll',n=300)
#' estpackG <- pack(x=x,greedy=TRUE)  ## estimate the packing number by greedy method
#' print(estpackG)
#' estpackC <- pack(x=x,greedy=FALSE) ## estimate the packing number by cluttering
#' print(estpackC)
#' @export
pack <- function(x, k1 = NULL, k2 = NULL, greedy = TRUE, eps = 0.01, DM = FALSE) {
    if (is.null(x)) {
        stop("data or distance matrix x is missing")
    }
    
    if (DM == FALSE) {
        Dmat <- as.matrix(stats::dist(x))
    } else {
        Dmat <- as.matrix(x)
    }
    n <- dim(Dmat)[1]
    
    if (greedy) {
        if (is.null(k1) || is.null(k2)) {
            R <- as.numeric(stats::quantile(as.numeric(Dmat[lower.tri(Dmat)]), seq(0, 1, len = 10))[2:3])
        } else {
            R <- c(k1, k2)
        }
        ## a greedy algorithm
        Lhist <- NULL
        for (l in 1:100) {
            ## permute data randomly
            ord <- sample(1:n, n, replace = FALSE)
            Lk <- NULL
            for (k in 1:2) {
                r <- R[k]
                Cset <- sample(1:n, 1)
                for (i in 2:n) {
                  CovID <- NULL
                  for (j in Cset) {
                    CovID <- c(CovID, which(Dmat[j, ] < r))
                  }
                  CovID <- unique(CovID)
                  if (length(CovID) < n) {
                    Cset <- rbind(Cset, setdiff(1:n, CovID)[1])
                  } else {
                    break
                  }
                }
                Lk <- c(Lk, log(length(Cset)))
            }
            Lhist <- rbind(Lhist, Lk)
            Dpack <- as.numeric(-(Lk[2] - Lk[1])/(log(R[2]) - log(R[1])))
            if (l > 10 && Dpack * (1 - eps)/2 > 1.65 * sqrt(sum(apply(Lhist, 2, stats::var)))/sqrt(l * 
                (log(R[2]/R[1])))) {
                return(max(Dpack, 0))
            }
        }
        return(max(Dpack, 0))
        
    } else {
        ## clustering for estimating the packing number
        if (is.null(k1) || is.null(k2)) {
            R <- as.numeric(stats::quantile(as.numeric(Dmat[lower.tri(Dmat)]), seq(0, 1, len = 10))[c(5, 
                10)])
        } else {
            R <- c(k1, k2)
        }
        That <- stats::hclust(stats::as.dist(Dmat))
        CLS <- stats::cutree(That, h = That$height)
        ANT <- CLS - CLS
        for (i in 1:dim(CLS)[2]) {
            cls <- CLS[, i]
            ncls <- length(unique(cls))
            for (j in 1:ncls) {
                MaxDist <- max(as.matrix(Dmat)[which(cls == j), which(cls == j)])
                ANT[which(cls == j), i] <- MaxDist
            }
        }
        Lk <- length(unique(CLS[, which(colMeans(ANT > R[1]) == 1)[1]]))
        Lk <- c(Lk, length(unique(CLS[, which(colMeans(ANT > R[2]) == 1)[1]])))
        
        Lk <- log(Lk)
        Dpack <- as.numeric(-(Lk[2] - Lk[1])/(log(R[2]) - log(R[1])))
        return(max(Dpack, 0))
    }
}
