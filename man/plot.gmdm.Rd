\name{plot.gmdm}
\alias{plot.gmdm}
\title{S3 method for class `gmdm'}
\usage{
  \method{plot}{gmdm} (x, labels, colors, main,
    ylab = "Fraction", xlab = "Position",
    label.length.max = 8, label.line.max = 3, cex.text = 1,
    cex.tickmark = 0.75, if.plot.new = TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{gmdm}.}

  \item{labels}{a string vector of the same length as
  \code{x$data}, giving the names of the numeric vectors in
  \code{x$data}.}

  \item{colors}{the colors of the discrete distributions;
  the default is \emph{"Dark2" colors in ColorBrewer
  palettes} if not specified.}

  \item{main}{an overall title for the plot. See
  \code{help("title", package="graphics")}; the default
  title is used if not specified.}

  \item{ylab}{a title for the y axis. See
  \code{help("title", package="graphics")}.}

  \item{xlab}{a title for the x axis. See
  \code{help("title", package="graphics")}.}

  \item{label.length.max}{numeric, giving the maximum
  string width allowed in diagonal labels.}

  \item{label.line.max}{numeric, giving the maximum number
  of lines allowed in diagonal labels.}

  \item{cex.text}{a numerical value giving the amount by
  which plot text should be magnified relative to the
  default.}

  \item{cex.tickmark}{a numerical value giving the amount
  by which tickmarks should be magnified relative to the
  default.}

  \item{if.plot.new}{logical, indicating whether to start a
  new plot device.}

  \item{...}{arguments to be passed to methods, see
  \code{gmdp}.}
}
\description{
  S3 method for class \code{gmdm}
}
\details{
  S3 method for class \code{gmdm}
}
\examples{
## ------------------------------------------------------------------------
## Example1: CAGE
## ------------------------------------------------------------------------
require("GMD") # load library
data(cage)     # load data

## construct a distance matrix and visualize it
short.labels <- gsub("(.+) \\\\(.+","\\\\1",names(cage)) # get short labels
x <- gmdm(cage[1:6],labels=short.labels[1:6])
plot(x)


## ------------------------------------------------------------------------
## Example2: ChIP-seq
## ------------------------------------------------------------------------
data(chipseq_mES)   # load data
data(chipseq_hCD4T) # load data

## pairwise distance and alignment based on GMD metric
plot(gmdm(chipseq_mES,sliding=FALSE))

## clustering on spatial distributions of histone modifications
x <- gmdm(chipseq_hCD4T,sliding=FALSE,resolution=10)
heatmap.3(x,revC=TRUE)
}
\references{
  See \code{help(GMD)}
}
\seealso{
  \code{\link{gmdm}}, \code{\link{gmdp}}
}
\keyword{hplot}
\keyword{methods}

