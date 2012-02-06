\name{plot.gmdp}
\alias{plot.gmdp}
\title{Plot function for class gmdp}
\usage{
  \method{plot}{gmdp} (x, labels = NULL, colors = NULL,
    main, ylab = "Fraction", xlab = "Position",
    xlim = NULL, if.text.gmd = TRUE, if.text.gap = TRUE,
    ...)
}
\arguments{
  \item{x}{an object of class \code{gmdp}.}

  \item{labels}{a string vector of the same length of
  \code{x$labels}, giving the names of the numeric vectors
  in \code{x}.}

  \item{colors}{the colors of the discrete distributions.
  See \code{help("plot.mhist", package="GMD")}.}

  \item{main}{an overall title for the plot.}

  \item{ylab}{a title for the y axis. See
  \code{help("plot.mhist", package="GMD")}.}

  \item{xlab}{a title for the x axis. See
  \code{help("plot.mhist", package="GMD")}.}

  \item{xlim}{numeric vectors of length 2, giving the x
  coordinates ranges.}

  \item{if.text.gmd}{logical, indicating whether
  \emph{GM-Distance} is reported in the subtitle.}

  \item{if.text.gap}{logical, indicating whether \emph{gap}
  is reported in the subtitle.}

  \item{...}{arguments to be passed to methods. See
  \code{help("plot.mhist", package="GMD")}.}
}
\description{
  Plot Function for Class \code{gmdp}
}
\details{
  Plot Function for Class \code{gmdp}
}
\examples{
require("GMD") # load library
data(cage)     # load data

## measure pairwise distance
x <- gmdp(cage[["Pfkfb3 (T02R00AEC2D8)"]],cage[["Csf1 (T03R0672174D)"]])
print(x)                     # print a brief version by default
print(x, mode="full")  # print a full version by default

## show alignment
plot(x,labels=c("Pfkfb3","Csf1"),beside=FALSE)

## show another alignment
plot(gmdp(cage[["Hig1 (T09R0743763C)"]],cage[["Cd72 (T04R028B8BC9)"]]),
     labels=c("Hig1 (T09R0743763C)","Cd72 (T04R028B8BC9)"),
     beside=FALSE)
}
\references{
  See \code{help(GMD)}
}
\seealso{
  \code{\link{gmdp}}
}
\keyword{hplot}
\keyword{methods}

