\name{gmdm}
\title{Generalized Minimum Distance Matrix}
\usage{gmdm(data,labels,pseudocount=0,sliding=TRUE,resolution=1)

## S3 method for class `gmdm'
print.gmdm(x, ...)

## convert a `gmdm' object into a `dist' object
gmdm2dist(m, diag=FALSE, upper=FALSE)

## compute GMDM and convert into a `dist' object
gmdm_dist(data, diag=FALSE, upper=FALSE, ...)}
\description{Computing Generalized Minimum Distance Matrix}
\details{Computing Generalized Minimum Distance Matrix}
\alias{gmdm}
\alias{print.gmdm}
\alias{gmdm2dist}
\alias{gmdm_dist}
\value{\code{gmdm} returns an object of class \code{gmdm}, a list with components

labels: a string vector, giving the names of distributions

data.ori: a list of numeric vectors, giving the original input

data: a list of numeric vectors, giving the normalized version of the original input

dm: a numeric numeric, the pairwise distance matrix of \emph{GM-Distances}

gap.pair: a numeric matrix, giving the gap pair of each alignment per row:
i.e. relative shifts between distributions of the optimal hit

sliding: logical, indicating whether sliding is performed

pseudocount: a numeric value that is allocated at each position in addition to original values}
\references{See \code{citation("GMD")}}
\seealso{\code{\link{plot.gmdm}}, \code{\link{gmdp}}}
\keyword{classes}
\arguments{\item{data}{a list of numeric vectors, a numeric matrix or data.frame}
\item{x}{a \code{gmdm} object.}
\item{m}{a \code{gmdm} object.}
\item{labels}{a character vector of the same length of x, giving the names of the numeric vectors.}
\item{pseudocount}{a numeric value to be allocated for each position to reduce bias;
by default \code{pseudocount = 0}.}
\item{sliding}{logical, indicating whether sliding is allowed or not for an optimal solution; 
by default \code{sliding = TRUE}.}
\item{resolution}{relative resolution, numeric (>=1), changing the size of the bin by
multiplying the value. A larger value (lower resolution) is more computational efficiet but
missing details.}
\item{diag}{logical value indicating whether the diagonal of the distance matrix should be 
printed by \code{print.dist}.}
\item{upper}{logical value indicating whether the upper triangle of the distance matrix should be
printed by \code{print.dist}.}
\item{...}{arguments to be passed to method}}
