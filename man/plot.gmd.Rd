\name{plot.gmd}
\alias{plot.gmd}
\title{Plot Function for Class gmd}
\usage{plot.gmd(x, labels=x$labels, colors=c("red", "blue"),
    plot.method=c("separate", "overlay"), plot.type="h", main,
    ylab="Fraction", xlab="Position", ylim, xlim, font.type=1,
    font.family=c("sans", "serif", "mono"), cex.lab=1.2,
    cex.tickmark=1, cex.legend=1.5, lwd.line=1, if.plot.new=TRUE,
    if.text.gmd=TRUE, if.text.gap=FALSE, if.plot.gap=TRUE,
    if.plot.legned=TRUE, x.jitter=ifelse(plot.method == "overlay",
    1/1000, 0), ...)}
\description{Plot Function for Class \code{gmd}}
\details{Plot Function for Class \code{gmd}}
\references{See \code{help(GMD)}}
\seealso{\code{\link{gmd}}}
\keyword{methods}
\keyword{hplot}
\arguments{\item{x}{an object of class \code{gmd}.}
\item{labels}{a string vector of the same length of \code{x$labels},
giving the names of the numeric vectors in \code{x}.}
\item{colors}{the colors of the discrete distributions; by default they are in \code{"red"} and \code{"blue"}.}
\item{plot.method}{the plot method. This can be specified as a string:
\code{"separate"}: means separated subplots [\emph{default}];
\code{"overlay"}: means overlaid subplots.}
\item{plot.type}{the plot type. See \code{type} in \code{plot} for possible values; the default \code{plot.type = "h"}, giving \sQuote{\bold{h}istogram} like vertical lines.}
\item{main}{an overall title for the plot. See \code{help("title", package="graphics")}.}
\item{ylab}{a title for the y axis. See \code{help("title", package="graphics")}.}
\item{xlab}{a title for the x axis. See \code{help("title", package="graphics")}.}
\item{ylim}{range of y values, as in \code{help("plot", package="graphics")}.}
\item{xlim}{range of x values, as in \code{help("plot", package="graphics")}.}
\item{font.type}{the name of a font type for drawing text. See \code{font} in \code{par}; the default \code{font.type = 1}, corresponding to plain text.}
\item{font.family}{the name of a font family for drawing text. See \code{family} in \code{par}; the default \code{font.family = "sans"}, corresponding to san serif typeface.}
\item{cex.lab}{a numerical value giving the amount by which \code{xlab} and \code{ylab} should be magnified relative to the default.}
\item{cex.tickmark}{a numerical value giving the amount by which tickmarks should be magnified relative to the default.}
\item{cex.legend}{a numerical value giving the amount by which legends should be magnified relative to the default.}
\item{lwd.line}{the line width, a \emph{positive} number, defaulting to \code{1}.}
\item{if.plot.new}{logical, indicating whether to start a new plot device.}
\item{if.text.gmd}{logical, indicating whether \emph{GM-Distance} is reported in the subtitle.}
\item{if.text.gap}{logical, indicating whether \emph{gap} is reported in the subtitle.}
\item{if.plot.gap}{logical, indicating whether \emph{gap} is plotted.}
\item{if.plot.legned}{logical, indicating whether \emph{legend} is plotted.}
\item{x.jitter}{numeric, indicating how \emph{jitter} should be added to distinguish subplots; by default \code{x.jitter=ifelse(plot.method=="overlay",1/1000,0)} giving how jitter should be adjusted according to the \emph{x-axis range}.}
\item{...}{arguments to be passed to methods, such as graphical parameters (see \code{par}).}}
\examples{require(GMD)
data(cage)
\dontrun{plot(gmd(cage[[1]],cage[[2]],labels=names(cage)[c(1,2)],
pseudocount=1, sliding=TRUE))}
plot(gmd(cage[[1]],cage[[3]],labels=names(cage)[c(1,3)],
pseudocount=1, sliding=TRUE))
plot(gmd(cage[[1]],cage[[3]],labels=names(cage)[c(1,3)],
pseudocount=1, sliding=TRUE), plot.method="overlay")}
