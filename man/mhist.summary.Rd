\name{mhist.summary}
\alias{mhist.summary}
\alias{plot.mhist.summary}
\title{Bin-wise summary of histograms}
\usage{
  mhist.summary(h, ...)

  \method{plot}{mhist.summary}(x,bins,plot.ci=TRUE,col=NULL,
    ci.color="orchid1",tcl=-0.25,omi=c(0.5,0.5,1.0,0.25),mar=c(3,3,3,1),
    mgp=c(2,0.5,0),if.plot.new=TRUE,...)
}
\arguments{
  \item{h}{a \code{"mhist"} object as produced by
  \code{as.mhist}}

  \item{x}{a \code{mhist.summary} object as produced by
  \code{mhist.summary}}

  \item{bins}{character vecter, the bin labels; if
  non-specific, bins are numbered/labeled starting with
  one.}

  \item{plot.ci}{logical, indicating whether plot error
  bars that represent the 0.50 confidence interval (CI)}

  \item{col}{color of the histogram}

  \item{ci.color}{color of the error bars}

  \item{tcl}{the length of tick marks as a fraction of the
  height of a line of text. See option \code{tcl}
  in\code{help("par", package="graphics")}.}

  \item{omi}{a vector of the form 'c(bottom, left, top,
  right)' giving the size of the outer margins in inches.
  See option \code{omi} in\code{help("par",
  package="graphics")}.}

  \item{mar}{a numerical vector of the form \code{c(bottom,
  left, top, right)} which gives the number of lines of
  margin to be specified on the four sides of the plot. See
  option \code{mar} in\code{help("par",
  package="graphics")}.}

  \item{mgp}{the margin line (in 'mex' units) for the axis
  title, axis labels and axis line.}

  \item{if.plot.new}{logical, whether starting a new device
  or not.}

  \item{...}{arguments to be passed to method
  \code{plot.mhist.summary}.\cr See \code{help("barplot2",
  package="gplots")}.}
}
\value{
  \code{mhist.summary} returns a \code{mhist.summary}
  object
}
\description{
  Bin-wise summary of a \code{mhist} object of histograms
}
\details{
  Bin-wise summary of a \code{mhist} object of histograms
}
\seealso{
  \code{\link{mhist}} \code{\link{plot.mhist}}
  \code{\link{plot.gmdp}} \code{\link{plot.gmdm}}
}

