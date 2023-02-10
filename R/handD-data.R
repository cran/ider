#' Hand ratation data
#'
#' Data from a QTL experiment on gravitropism in
#' Arabidopsis, with data on 162 recombinant inbred lines (Ler x
#' Cvi). The outcome is the root tip angle (in degrees) at two-minute
#' increments over eight hours.
#'
#' @docType data
#'
#' @usage data(handD)
#'
#' @format An object of class \code{'dist'}.
#'
#' @keywords datasets
#'
#' @references  E. Levina and P. J. Bickel. Maximum likelihood estimation of intrinsic dimension. Advances in Neural Information Processing Systems 17, 2005.
#' @references B. Kegl. Intrinsic dimension estimation using packing numbers. Advances in Neural Information Processing Systems 15, 2002.
#' @references H. Hino, J. Fujiki, S. Akaho, and N. Murata, 'Local Intrinsic Dimension Estimation by Generalized Linear Modeling', Neural Computation, 2017
#'
#'
#' @examples
#' data(handD)
#' estmle <- lbmle(x=handD,DM=TRUE,k1=5,k2=10)
#' print(estmle)
"handD"
