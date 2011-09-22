\name{gmdm}
\alias{gmdm}
\title{Generalized Minimum Distance Matrix}
\usage{gmdm(x, labels=names(x), pseudocount=0, sliding=TRUE)}
\description{Generalized Minimum Distance Matrix}
\details{Generalized Minimum Distance Matrix}
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
\seealso{\code{\link{plot.gmdm}}, \code{\link{gmd}}}
\keyword{classes}
\arguments{\item{x}{a list of numeric vectors}
\item{labels}{a character vector of the same length of x, giving the names of the numeric vectors.}
\item{pseudocount}{a numeric value to be allocated for each position to reduce bias;
by default \code{pseudocount = 0}.}
\item{sliding}{logical, indicating whether sliding is allowed or not for an optimal solution; 
by default \code{sliding = TRUE}.}}
\examples{require(GMD)
data(cage)
x <- gmdm(cage)
print(x$labels)
print(x$dm)

\dontrun{data(cagel)
x <- gmdm(cagel)
head(x$labels)
head(x$dm)}}
