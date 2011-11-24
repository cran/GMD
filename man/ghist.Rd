\name{ghist}
\title{Generalized Histogram Computation}
\usage{ghist(data,n,breaks=if (!invalid(n)) NULL else "Sturges",
bins=NULL,digits=1)

gbreaks(data, n)

is.ghist(x)

as.ghist(x,bins)

is.mhist(x)

as.mhist(x,bins)

mhist2matrix(h)}
\description{Generalized Histogram Computation with classes to contain a single histogram or multiple histograms}
\details{\code{ghist} generates a single histogram.

\code{gbreaks} generate bin boundaries for a histogram.

\code{is.ghist} returns \code{TRUE} if \code{x} is an object of code{ghist} and
\code{FALSE} otherwise.

\code{as.ghist} is a generic function.  The method for numeric vectors will return a
\code{ghist} object.

\code{is.mhist} returns \code{TRUE} if \code{x} is an object of code{mhist} and
\code{FALSE} otherwise.

\code{as.mhist} is a generic function.  The method is for numeric list, matrices or data frames and will return a
\code{mhist} object.

\code{mhist2matrix} convert a \code{mhist} object into a numeric matrix, filling observations by row.}
\alias{ghist}
\alias{gbreaks}
\alias{is.ghist}
\alias{as.ghist}
\alias{mhist}
\alias{is.mhist}
\alias{as.mhist}
\alias{mhist2matrix}
\seealso{\code{\link{plot.mhist}} \code{\link{mhist.summary}} \code{\link{plot.mhist.summary}}}
\arguments{\item{data}{a vector of values for which the histogram is desired.}
\item{n}{a single number giving the number of bins for the histogram.}
\item{breaks}{a vector giving the breakpoints between histogram bins, or
a character string naming an algorithm to compute the number of bins, or
a function to compute the number of bins (see \code{help("dist", package="graphics")}).}
\item{bins}{character vecter, the bin labels.}
\item{digits}{integer, the number of digits to round for breaks.}
\item{x}{an R object.}
\item{h}{an object of class \code{mhist}}}
\examples{## load library
require("GMD")

## create two normally-distributed samples
## with unequal means and unequal variances
set.seed(2012)
v1 <- rnorm(1000,mean=-5, sd=10)
v2 <- rnorm(1000,mean=10, sd=5)

## create common bins
n <- 20 # desired number of bins
breaks <- gbreaks(c(v1,v2),n) # bin boundaries
x <-
list(ghist(v1,breaks=breaks,digits=0),
ghist(v2,breaks=breaks,digits=0))
mhist.obj <- as.mhist(x)}
