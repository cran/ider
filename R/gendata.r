#' Data generator for intrinsic dimension estimation.
#'
#' \code{gendata} generates various artificial datasets for intrinsic dimension estimation experiments.
#' 
#' This function generates various artificial datasets often used in 
#' manifold learning and dimension estimation researches.
#' For some datasets, complexity of the shape is controlled by the parameter \code{curv}.
#' The parameters \code{noise} and \code{outlier} are used for adding noise and/or 
#' outliers for the dataset.
#'
#'
#' @param DataName Name of dataset, one of the following: 
#' \itemize{
#'  \item SwissRoll: SwissRoll data, 2D manifold in 3D space.
#'  \item NDSwissRoll: Non-deformable SwissRoll data, 2D manifold in 3D space.
#'  \item Moebius: Moebius strip, 2D manifold in 3D space.
#'  \item SphericalShell: Spherical Shell, (p-1)-dimensional manifold in p-dimensional space.
#'  \item Sinusoidal: Sinusoidal data, 1D manifold in 3D space.
#'  \item Spiral: Spiral-shaped 1D manifold in 2D space.
#'  \item Cylinder: Cylinder-shaped 2D manifold in 3D space.
#'  \item SShape: S-shaped 2D manifold in 3D space.
#'  \item ldbl: LDB(line - disc - filled ball - line), embedded in 3D space (original dataset).
#' }
#' @param n number of data points to be generated.
#' @param p ambient dimension of the dataset.
#' @param noise parameter to control noise level in the dataset. In many cases,
#' it is used for \code{sd} of \code{rnorm} used inside the function.
#' @param ol percentage of outliers, i.e., n * ol outliers are added to the generated dataset.
#' @param curv a parameter to control the complexity of the embedded manifold.
#' @param seed random number seed.
#' @param sorted logical. If \code{TRUE}, the index of the generated dataset is sorted
#' with respect to x-axis for the ease of visualization.
#' @return Data matrix. For \code{ldbl} dataset, it outputs a list composed of
#'  \code{x}: data matrix and \code{tDim}: true intrinsic dimension for each point.
#' @author Hideitsu Hino \email{hideitsu.hino@@gmail.com}
#' @examples 
#' ## global intrinsic dimension estimate
#' x <- gendata(DataName='SwissRoll')
#' estmle <- lbmle(x=x,k1=3,k2=5)
#' print(estmle)
#' 
#' ## local intrinsic dimension estimate
#' tmp <- gendata(DataName='ldbl',n=1000)
#' x <- tmp$x
#' estmada <- mada(x=x,local=TRUE)
#' head(estmada)  ## estimated local intrinsic dimensions
#' head(tmp$tDim) ## true local intrinsic dimensions
#' @export
gendata <- function(DataName = "SwissRoll", n = 300, p = NULL, noise = NULL, ol = NULL, curv = 1, 
    seed = 123, sorted = FALSE) {
    set.seed(seed)
    
    if (DataName == "Moebius") {
        p <- 3
        q <- 2
        u <- stats::runif(n, -1, 1)
        v <- stats::runif(n, 0, 2 * pi)
        if (sorted) {
            v <- sort(v)
        }
        ret <- cbind((1 + 0.5 * u * cos(0.5 * v * curv)) * cos(v), (1 + 0.5 * u * cos(0.5 * 
            v * curv)) * sin(v), 0.5 * u * sin(0.5 * v * curv))
        if (!is.null(noise)) {
            ret <- ret + matrix(stats::rnorm(n * p, m = 0, sd = noise), n, p)
        }
        if (!is.null(ol)) {
            no <- ceiling(n * ol)
            tmp <- matrix(stats::runif(no * p, range(ret)[1], range(ret)[2]), no, p)
            ret[sample(dim(ret)[1], no), ] <- tmp
        }
        return(ret)
    } else if (DataName == "Spiral") {
        p <- 2
        q <- 1
        x <- stats::runif(n, 0, 1)
        if (sorted) {
            x <- sort(x)
        }
        ret <- cbind(sqrt(x) * cos(10 * pi * sqrt(x)), sqrt(x) * sin(10 * pi * sqrt(x)))
        ret <- ret + matrix(stats::rnorm(2 * n, m = 0, sd = 0.01), ncol = 2)
        
        if (!is.null(noise)) {
            ret <- ret + matrix(stats::rnorm(n * p, m = 0, sd = noise), n, p)
        }
        if (!is.null(ol)) {
            no <- ceiling(n * ol)
            tmp <- matrix(stats::runif(no * p, range(ret)[1], range(ret)[2]), no, p)
            ret[sample(dim(ret)[1], no), ] <- tmp
        }
        return(ret)
    } else if (DataName == "SphericalShell") {
        if (p <= 1) {
            stop("Too small embedding dimension")
        }
        q <- p - 1
        x1 <- stats::runif(n * p, -1, 1)
        tmp <- matrix(x1, n, p)
        ret <- apply(tmp, 2, FUN = function(x) {
            x/sqrt(rowSums(tmp^2))
        })
        
        if (sorted) {
            ret <- ret[order(ret[, 3]), ]
        }
        
        if (!is.null(noise)) {
            ret <- ret + matrix(stats::rnorm(n * p, m = 0, sd = noise), n, p)
        }
        if (!is.null(ol)) {
            no <- ceiling(n * ol)
            tmp <- matrix(stats::runif(no * p, range(ret)[1], range(ret)[2]), no, p)
            ret[sample(dim(ret)[1], no), ] <- tmp
        }
        return(ret)
        
    } else if (DataName == "Sinusoidal") {
        p <- 3
        q <- 1
        u <- stats::runif(n, 0, 2 * pi)
        ret <- cbind(sin(u), cos(u), 0.1 * sin(curv * u))
        if (!is.null(noise)) {
            ret <- ret + matrix(stats::rnorm(n * p, m = 0, sd = noise), n, p)
        }
        if (!is.null(ol)) {
            no <- ceiling(n * ol)
            tmp <- matrix(stats::runif(no * p, range(ret)[1], range(ret)[2]), no, p)
            ret[sample(dim(ret)[1], no), ] <- tmp
        }
        return(ret)
        
    } else if (DataName == "Trefoil") {
        p <- 3
        q <- 1
        u <- stats::runif(n, 0, 2 * pi)
        if (sorted) {
            u <- sort(u)
        }
        ret <- cbind(41 * cos(u) - 18 * sin(u) - 83 * cos(u) - 83 * sin(2 * u) - 11 * sin(3 * 
            u) + 27 * sin(3 * u), 36 * cos(u) + 27 * sin(u) - 113 * cos(2 * u) + 30 * sin(2 * 
            u) + 11 * cos(3 * u) - 27 * sin(3 * u), 45 * sin(u) - 30 * cos(2 * u) + 113 * sin(2 * 
            u) - 11 * cos(3 * u) + 27 * sin(3 * u))
        
        if (!is.null(noise)) {
            ret <- ret + matrix(stats::rnorm(n * p, m = 0, sd = noise * mean(abs(range(ret)))), 
                n, p)
        }
        if (!is.null(ol)) {
            no <- ceiling(n * ol)
            tmp <- matrix(stats::runif(no * p, range(ret)[1], range(ret)[2]), no, p)
            ret[sample(dim(ret)[1], no), ] <- tmp
        }
        return(ret)
        
    } else if (DataName == "Torus") {
        p <- 3
        q <- 2
        x1 <- stats::runif(n, 0, 2 * pi)
        x2 <- stats::runif(n, 0, 2 * pi)
        if (sorted) {
            x1 <- sort(x1)
        }
        ret <- cbind((2 + cos(x1)) * cos(x2), (2 + cos(x1)) * sin(x2), sin(x1))
        
        if (!is.null(noise)) {
            ret <- ret + matrix(stats::rnorm(n * p, m = 0, sd = noise * mean(abs(range(ret)))), 
                n, p)
        }
        if (!is.null(ol)) {
            no <- ceiling(n * ol)
            tmp <- matrix(stats::runif(no * p, range(ret)[1], range(ret)[2]), no, p)
            ret[sample(dim(ret)[1], no), ] <- tmp
        }
        
        return(ret)
        
    } else if (DataName == "SShape") {
        p <- 3
        q <- 2
        x <- seq(-1, 1, len = n)
        y <- stats::runif(n, -1, 1)
        ret <- cbind(x, y, sin(pi * x * curv))
        if (!is.null(noise)) {
            ret <- ret + matrix(stats::rnorm(n * p, m = 0, sd = noise), n, p)
        }
        
        if (!is.null(ol)) {
            no <- ceiling(n * ol)
            tmp <- matrix(stats::runif(no * p, range(ret)[1], range(ret)[2]), no, p)
            ret[sample(dim(ret)[1], no), ] <- tmp
            
        }
        return(ret)
        
    } else if (DataName == "Cylinder") {
        p <- 3
        q <- 2
        x <- sqrt(2 + 2 * seq(-1, 1 - 2/n, 2/n))
        y <- 2 * stats::runif(n, -1, 1)
        ret <- cbind(cos(2 * pi * x * curv), y, sin(2 * pi * x * curv))
        if (!is.null(noise)) {
            ret <- ret + matrix(stats::rnorm(n * p, m = 0, sd = noise), n, p)
        }
        
        if (!is.null(ol)) {
            no <- ceiling(n * ol)
            tmp <- matrix(stats::runif(no * p, range(ret)[1], range(ret)[2]), no, p)
            ret[sample(dim(ret)[1], no), ] <- tmp
        }
        return(ret)
        
    } else if (DataName == "SwissRoll") {
        p <- 3
        q <- 2
        x1 <- stats::runif(n, -1, 1)
        if (sorted) {
            x1 <- sort(x1)
        }
        y <- stats::runif(n, -1, 1)
        ret <- cbind(sqrt(2 + 2 * x1) * cos(2 * pi * sqrt(2 + 2 * x1)), sqrt(2 + 2 * x1) * 
            sin(2 * pi * sqrt(2 + 2 * x1)), 2 * y)
        if (!is.null(noise)) {
            ret <- ret + matrix(stats::rnorm(n * p, m = 0, sd = noise), n, p)
        }
        
        if (!is.null(ol)) {
            no <- ceiling(n * ol)
            tmp <- matrix(stats::runif(no * p, range(ret)[1], range(ret)[2]), no, p)
            ret[sample(dim(ret)[1], no), ] <- tmp
        }
        return(ret)
        
    } else if (DataName == "NDSwissRoll") {
        p <- 3
        q <- 2
        x1 <- stats::runif(n, -1, 1)
        if (sorted) {
            x1 <- sort(x1)
        }
        y <- stats::runif(n, -1, 1)
        ret <- cbind((1 + y^(2)) * sqrt(1 + x1) * cos(2 * pi * sqrt(1 + x1)), (1 + y^(2)) * 
            sqrt(1 + x1) * sin(2 * pi * sqrt(1 + x1)), 2 * y)
        if (!is.null(noise)) {
            ret <- ret + matrix(stats::rnorm(n * p, m = 0, sd = noise), n, p)
        }
        
        if (!is.null(ol)) {
            no <- ceiling(n * ol)
            tmp <- matrix(stats::runif(no * p, range(ret)[1], range(ret)[2]), no, p)
            ret[sample(dim(ret)[1], no), ] <- tmp
        }
        return(ret)
        
    } else if (DataName == "ldbl") {
        ## line - disc - filled ball - line, again
        
        ## line in 3D
        line <- cbind(rep(0, 5 * n), rep(0, 5 * n), stats::runif(5 * n, -0.5, 0))
        
        ## disc in 3D
        disc <- cbind(matrix(stats::runif(26 * n, -1, 1), ncol = 2), rep(0, n))
        disc <- disc[-which(sqrt(rowSums(disc^2)) > 1), ]
        disc <- disc[, c(1, 3, 2)]
        disc[, 3] <- disc[, 3] - min(disc[, 3]) + max(line[, 3])
        
        ## filled ball in 3D
        fb <- matrix(stats::runif(30 * n * 10, -0.5, 0.5), ncol = 3)
        
        rmID <- which(sqrt(rowSums(fb^2)) > 0.5)
        if (length(rmID) > 0) {
            fb <- fb[-which(sqrt(rowSums(fb^2)) > 0.5), ]
        }
        fb <- cbind(fb[, 1:2], fb[, 3] + 0.5)
        
        fb[, 3] <- fb[, 3] - min(fb[, 3]) + max(disc[, 3])
        
        if (sorted) {
            fb <- fb[order(fb[, 3]), ]
        }
        
        ## line in 3D
        line2 <- cbind(rep(0, 5 * n), rep(0, 5 * n), stats::runif(5 * n, -0.5, 0))
        line2[, 3] <- line2[, 3] - min(line2[, 3]) + max(fb[, 3])
        
        lineID <- rep(1, dim(line)[1])
        discID <- rep(2, dim(disc)[1])
        fbID <- rep(3, dim(fb)[1])
        line2ID <- rep(1, dim(line2)[1])
        x <- rbind(line, disc, fb, line2)
        useID <- sort(sample(1:dim(x)[1], n))
        x <- x[useID, ]
        
        return(list(x = x, tDim = c(lineID, discID, fbID, line2ID)[useID]))
        
    } else {
        stop("Invalid Data Name")
    }
    
}
