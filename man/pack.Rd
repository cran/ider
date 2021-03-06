% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pack.r
\name{pack}
\alias{pack}
\title{Intrinsic Dimension Estimation Using Packing Numbers.}
\usage{
pack(x, k1 = NULL, k2 = NULL, greedy = TRUE, eps = 0.01, DM = FALSE)
}
\arguments{
\item{x}{data matrix or distance matrix given by as.matrix(dist(x)).}

\item{k1}{first radius parameter. If one of \code{k1} or \code{k2} is \code{NULL},
then both are automatically determined from the input data.}

\item{k2}{second radius parameter.}

\item{greedy}{logical. If \code{TRUE}, then a greedy algorithm is used for 
estimating the packing number. If \code{FALSE}, then a hierarchical clustering
algorithm is used instead.}

\item{eps}{accuracy parameter.}

\item{DM}{whether \code{'x'} is distance matrix or not. logical.}
}
\value{
Estimated global intrinsic dimension.
}
\description{
\code{pack} estimates intrinsic dimension of given dataset based on the packing number.
}
\details{
A variant of fractal dimension called the capacity dimension is considered.
The capacity dimension is defined by using the notion of covering number, which is
hard to calculate in general. In this function, the packing number of the data
 space is used as the surrogate of the covering number.
  The packing number is estimated by greedy manner or by hierarchical clustering.
}
\examples{
x <- gendata(DataName='SwissRoll',n=300)
estpackG <- pack(x=x,greedy=TRUE)  ## estimate the packing number by greedy method
print(estpackG)
estpackC <- pack(x=x,greedy=FALSE) ## estimate the packing number by cluttering
print(estpackC)
}
\references{
B. Kegl. Intrinsic dimension estimation using packing numbers.
 Advances in Neural Information Processing Systems 15, 2002.

B. Eriksson and M. Crovella. Estimating intrinsic dimension via clustering.
IEEE Statistical Signal Processing Workshop, 2012.
}
\author{
Hideitsu Hino \email{hideitsu.hino@gmail.com}
}
