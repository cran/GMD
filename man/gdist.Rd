\name{gdist}
\title{Generalized Distance Matrix Computation}
\usage{gdist(x,method="euclidean",MoreArgs=NULL,diag=FALSE,upper=FALSE)

is.dist(d)}
\description{\code{gdist} computes and returns the distance matrix computed by using user-defined distance measure.}
\details{\code{is.dist} tests if its argument is a `dist' object.

The distance (or dissimilarity) function (\code{FUN}) can be any distance measure applied to \code{x}.
For instance, \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"},\code{"canberra"},
\code{"binary"}, \code{"minkowski"}, "correlation.of.variables", "correlation.of.observations" or
\code{gmdm}. "correlation.of.variables" computes the correlation distance of
the variables (the columns); all the other compute the distances between
the observations (the rows) of a data matrix.}
\alias{gdist}
\alias{is.dist}
\value{\code{gdist} returns an object of `dist'.\cr}
\arguments{\item{x}{a numeric matrix, data frame or `dist' object.}
\item{method}{the distance measure to be used. This can either be one of
the methods used in \code{dist} (see \code{help("dist", package="stats")})
or \code{"correlation"}, \code{"correlation.of.observations"} and
\code{"correlation.of.variables"}. In addition, user-defined distance measure
are also allowed, which returns a \emph{dist} object and should at least
have attributes \emph{"Size"} and \emph{"Labels"}.}
\item{MoreArgs}{a list of other arguments to be passed to \code{gdist}.}
\item{diag}{logical value indicating whether the diagonal of the distance matrix should be 
printed by \code{print.dist}.}
\item{upper}{logical value indicating whether the upper triangle of the distance matrix should be
printed by \code{print.dist}.}
\item{d}{an R object.}}
\examples{## load library
require("GMD")
require(cluster)

## compute distance using Euclidean metric (default)
data(ruspini)
x <- gdist(ruspini)

## see a dendrogram result by hierarchical clustering
dev.new(width=12, height=6)
plot(hclust(x),
main="Cluster Dendrogram of Ruspini data",
xlab="Observations")

## convert to a distance matrix
m <- as.matrix(x)

## convert from a distance matrix
d <- as.dist(m)
stopifnot(d == x)

## Use correlations between variables "as distance"
data(USJudgeRatings)
dd <- gdist(x=USJudgeRatings,method="correlation.of.variables")
dev.new(width=12, height=6)
plot(hclust(dd),
main="Cluster Dendrogram of USJudgeRatings data",
xlab="Variables")}
