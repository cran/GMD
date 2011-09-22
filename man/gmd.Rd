\name{gmd}
\alias{gmd}
\title{Generalized Minimum Distance (GMD)}
\usage{gmd(v1, v2, labels=c("v1", "v2"), pseudocount=0, sliding=TRUE)}
\description{Generalized Minimum Distance}
\details{Generalized Minimum Distance}
\value{\code{gmd} returns an object of class \code{gmd}, a list with components

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
\seealso{\code{\link{print.gmd}}, \code{\link{summary.gmd}}, \code{\link{plot.gmd}}, \code{\link{gmdm}}}
\keyword{classes}
\arguments{\item{v1}{a numeric vector, giving positional counts as a discrete distribution.}
\item{v2}{a numeric vector, giving positional counts as a discrete distribution.}
\item{labels}{a string vector of length 2, giving the names of v1 and v2 respectively.}
\item{pseudocount}{a numeric value to be allocated for each position to reduce bias;
by default \code{pseudocount = 0}.}
\item{sliding}{logical, indicating whether sliding is allowed or not for an optimal solution; 
by default \code{sliding = TRUE}.}}
\examples{require(GMD)
gmd(c(4,1,1,0,0,0,3,1),c(2,1,1,0,0,0,3,3),sliding=FALSE)
x <- gmd(c(4,1,1,0,0,0,3,1), c(1,1,2,1,1,0,0,0,3,3,5,5),
pseudocount=1, sliding=TRUE)
print(x)
print(x, "full")}
