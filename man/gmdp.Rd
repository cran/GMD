\name{gmdp}
\title{Generalized Minimum Distance between a pair of distributions}
\usage{gmdp(v1, v2, labels=c("v1","v2"), pseudocount=0, sliding=TRUE,
resolution=1)

## S3 method for class `gmdp'
print.gmdp(x, print.mode=c("brief","detailed","full"),
digits=3, ...)

## S3 method for class `gmdp'
summary.gmdp(object, ...)}
\description{Generalized Minimum Distance between a pair of distributions}
\details{Generalized Minimum Distance between a pair of distributions}
\alias{gmdp}
\alias{print.gmdp}
\alias{summary.gmdp}
\value{\code{gmdp} returns an object of class \code{gmdp}, a numeric with an attribute of
\emph{meta} in a list with components: 

labels: a string vector, giving the names of distributions

v1.ori: a numeric vector, the first input distribution

v2.ori: a numeric vector, the second input distribution

v1: a numeric vector, the normalized version of the first input distribution

v2: a numeric vector, the normalized version of the second input distribution

distance: numeric, the \emph{GM-Distance} (\emph{GMD})

sliding: logical, indicating whether sliding is performed

pseudocount: a numeric value that is allocated at each position in addition to original values

gap.pair: a numeric matrix, giving one gap pair per row:
i.e. relative shifts between distributions of one optimal hit

n.hit: numeric, the number of (equally good) optimal hits}
\references{See \code{citation("GMD")}}
\seealso{\code{\link{print.gmdp}}, \code{\link{summary.gmdp}}, \code{\link{plot.gmdp}},
\code{\link{gmdm}}}
\keyword{classes}
\arguments{\item{v1}{a numeric vector, giving positional counts as a discrete distribution.}
\item{v2}{a numeric vector, giving positional counts as a discrete distribution.}
\item{labels}{a string vector of length 2, giving the names of v1 and v2 respectively.}
\item{pseudocount}{a numeric value to be allocated for each position to reduce bias;
by default \code{pseudocount = 0}.}
\item{sliding}{logical, indicating whether sliding is allowed or not for an optimal solution; 
by default \code{sliding = TRUE}.}
\item{resolution}{relative resolution, numeric (>=1), changing the size of the bin by
multiplying the value. A larger value (lower resolution) is more computational efficiet
but missing details.}
\item{x}{an object of class \code{gmdp}.}
\item{object}{an object of class \code{gmdp}.}
\item{print.mode}{a string of the following: \code{c("brief","detailed","full")},
indicating whether to print in \emph{full} mode (\emph{default}).}
\item{digits}{integer, indicating the number of decimal places to be printed.}
\item{...}{arguments to be passed to method.}}
\examples{require(GMD)
gmdp(c(4,1,1,0,0,0,3,1),c(2,1,1,0,0,0,3,3),sliding=FALSE)
x <- gmdp(c(4,1,1,0,0,0,3,1), c(1,1,2,1,1,0,0,0,3,3,5,5),
pseudocount=1, sliding=TRUE)
print(x)
print(x, "full")}
