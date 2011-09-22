\name{plot.gmdm}
\alias{plot.gmdm}
\title{Plot Function for Class gmdm}
\usage{plot.gmdm(x, labels=x$labels, colors, plot.type="h", main,
    ylab="Fraction", xlab="Position", label.length.max=8,
    label.line.max=3, cex.text=2, cex.tickmark=0.75, if.plot.new=TRUE,
    x.jitter=1/1000, ...)}
\description{Plot Function for Class \code{gmdm}}
\details{Plot Function for Class \code{gmdm}}
\references{See \code{help(GMD)}}
\seealso{\code{\link{gmdm}}, \code{\link{gmd}}}
\keyword{methods}
\keyword{hplot}
\arguments{\item{x}{an object of class \code{gmdm}.}
\item{labels}{a string vector of the same length as \code{x$data},
giving the names of the numeric vectors in \code{x$data}.}
\item{colors}{the colors of the discrete distributions; the default is \emph{"Dark2" colors in ColorBrewer palettes} if not specified.}
\item{plot.type}{the plot type. See \code{type} in \code{plot} for possible values; the default \code{plot.type = "h"}, giving \sQuote{\bold{h}istogram} like vertical lines.}
\item{main}{an overall title for the plot. See \code{help("title", package="graphics")}; the default title is used if not specified.}
\item{ylab}{a title for the y axis. See \code{help("title", package="graphics")}.}
\item{xlab}{a title for the x axis. See \code{help("title", package="graphics")}.}
\item{label.length.max}{numeric, giving the maximum string width allowed in diagonal labels.}
\item{label.line.max}{numeric, giving the maximum number of lines allowed in diagonal labels.}
\item{cex.text}{a numerical value giving the amount by which plot text should be magnified relative to the default.}
\item{cex.tickmark}{a numerical value giving the amount by which tickmarks should be magnified relative to the default.}
\item{if.plot.new}{logical, indicating whether to start a new plot device.}
\item{x.jitter}{numeric, indicating how \emph{jitter} should be added to distinguish subplots; by default \code{x.jitter=1/1000} indicating the jitter is adjusted to 1/1000 of the \emph{x-axis range}.}
\item{...}{arguments to be passed to methods, see \code{gmd}.}}
\examples{require(GMD)
data(cage)
plot(gmdm(cage))}
