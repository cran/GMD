\name{legend}
\alias{legend}
\title{Add Legends to Plots}
\usage{
  legend(x, y = NULL, legend, fill = NULL,
    col = par("col"), border = "black", lty, lwd, pch,
    angle = 45, density = NULL, bty = "o", bg = par("bg"),
    box.lwd = par("lwd"), box.lty = par("lty"),
    box.col = par("fg"), pt.bg = NA, cex = 1, pt.cex = cex,
    pt.lwd = lwd, xjust = 0, yjust = 1, x.intersp = 1,
    y.intersp = 1, adj = c(0, 0.5), text.width = NULL,
    text.col = par("col"), merge = do.lines && has.pch,
    trace = FALSE, plot = TRUE, ncol = 1, horiz = FALSE,
    title = NULL, inset = 0, xpd, title.col = text.col,
    title.adj = 0.5, seg.len = 2)
}
\arguments{
  \item{x}{the x coordinates to be used to position the
  legend.}

  \item{y}{the y coordinates to be used to position the
  legend. \code{x} and \code{y} can be specified by keyword
  or in any way which is accepted by
  \code{\link{xy.coords}}: See \sQuote{Details}.}

  \item{legend}{a character or \link{expression} vector. of
  length \eqn{\ge 1}{>= 1} to appear in the legend.  Other
  objects will be coerced by
  \code{\link{as.graphicsAnnot}}.}

  \item{fill}{if specified, this argument will cause boxes
  filled with the specified colors (or shaded in the
  specified colors) to appear beside the legend text.}

  \item{col}{the color of points or lines appearing in the
  legend.}

  \item{border}{the border color for the boxes (used only
  if \code{fill} is specified).}

  \item{lty}{the line types for lines appearing in the
  legend.}

  \item{lwd}{the line widths for lines appearing in the
  legend. One of \code{lty} and \code{lwd} \emph{must} be
  specified for line drawing.}

  \item{pch}{the plotting symbols appearing in the legend,
  either as vector of 1-character strings, or one (multi
  character) string.  \emph{Must} be specified for symbol
  drawing.}

  \item{angle}{angle of shading lines.}

  \item{density}{the density of shading lines, if numeric
  and positive. If \code{NULL} or negative or \code{NA}
  color filling is assumed.}

  \item{bty}{the type of box to be drawn around the legend.
  The allowed values are \code{"o"} (the default) and
  \code{"n"}.}

  \item{bg}{the background color for the legend box.  (Note
  that this is only used if \code{bty != "n"}.)}

  \item{box.lwd}{the line type for the legend box.}

  \item{box.lty}{the line width for the legend box.}

  \item{box.col}{the color for the legend box.}

  \item{pt.bg}{the background color for the
  \code{\link{points}}, corresponding to its argument
  \code{bg}.}

  \item{cex}{character expansion factor \bold{relative} to
  current \code{par("cex")}.}

  \item{pt.cex}{expansion factor(s) for the points.}

  \item{pt.lwd}{line width for the points, defaults to the
  one for lines, or if that is not set, to
  \code{par("lwd")}.}

  \item{xjust}{how the legend is to be justified relative
  to the legend x location.  A value of 0 means left
  justified, 0.5 means centered and 1 means right
  justified.}

  \item{yjust}{the same as \code{xjust} for the legend y
  location.}

  \item{x.intersp}{character interspacing factor for
  horizontal (x) spacing.}

  \item{y.intersp}{the same for vertical (y) line
  distances.}

  \item{adj}{numeric of length 1 or 2; the string
  adjustment for legend text.  Useful for y-adjustment when
  \code{labels} are \link{plotmath} expressions.}

  \item{text.width}{the width of the legend text in x
  (\code{"user"}) coordinates.  (Should be positive even
  for a reversed x axis.) Defaults to the proper value
  computed by \code{\link{strwidth}(legend)}.}

  \item{text.col}{the color used for the legend text.}

  \item{merge}{logical; if \code{TRUE}, merge points and
  lines but not filled boxes.  Defaults to \code{TRUE} if
  there are points and lines.}

  \item{trace}{logical; if \code{TRUE}, shows how
  \code{legend} does all its magical computations.}

  \item{plot}{logical.  If \code{FALSE}, nothing is plotted
  but the sizes are returned.}

  \item{ncol}{the number of columns in which to set the
  legend items (default is 1, a vertical legend).}

  \item{horiz}{logical; if \code{TRUE}, set the legend
  horizontally rather than vertically (specifying
  \code{horiz} overrides the \code{ncol} specification).}

  \item{title}{a character string or length-one expression
  giving a title to be placed at the top of the legend.
  Other objects will be coerced by
  \code{\link{as.graphicsAnnot}}.}

  \item{inset}{inset distance(s) from the margins as a
  fraction of the plot region when legend is placed by
  keyword.}

  \item{xpd}{if supplied, a value of the graphical
  parameter 'xpd' to be used while the legend is being
  drawn.}

  \item{title.col}{color for \code{title}.}

  \item{title.adj}{horizontal adjustment for \code{title}:
  see the help for \code{par("adj")}.}

  \item{seg.len}{the length of lines drawn to illustrate
  \code{lty} and/or \code{lwd} (in units of character
  widths).}
}
\description{
  This function can be used to add legends to plots.  Note
  that a call to the function \code{\link{locator}(1)} can
  be used in place of the \code{x} and \code{y} arguments.
}
\details{
  see \code{legend} in package:graphics for details; Note:
  Old versions of graphics:::legend do not have `border'
  option.
}

