\name{heatmap.3}
\alias{heatmap.3}
\title{Enhanced Heatmap Representation with Dendrogram and Partition}
\usage{
  heatmap.3(x, diss = inherits(x, "dist"), Rowv = TRUE,
    Colv = TRUE,
    dendrogram = c("both", "row", "column", "none"),
    dist.row, dist.col, dist.FUN = gdist,
    dist.FUN.MoreArgs = list(method = "euclidean"),
    hclust.row, hclust.col, hclust.FUN = hclust,
    hclust.FUN.MoreArgs = list(method = "ward"),
    scale = c("none", "row", "column"), na.rm = TRUE,
    cluster.by.row = FALSE, cluster.by.col = FALSE,
    kr = NA, kc = NA, row.clusters = NA, col.clusters = NA,
    revR = FALSE, revC = FALSE, add.expr, breaks, x.center,
    color.FUN = "bluered", sepList = list(NULL, NULL),
    sep.color = c("gray45", "gray45"), sep.lty = 1,
    sep.lwd = 2, cellnote, cex.note = 1, notecol = "cyan",
    na.color = par("bg"),
    trace = c("none", "column", "row", "both"),
    tracecol = "cyan", hline, vline, linecol = tracecol,
    labRow = TRUE, labCol = TRUE, margin.for.labRow,
    margin.for.labCol, ColIndividualColors,
    RowIndividualColors, cexRow, cexCol,
    labRow.by.group = FALSE, labCol.by.group = FALSE,
    key = TRUE, keysize = 1.5, mapsize = 9, mapratio = 4/3,
    sidesize = 3, cex.key.main = 0.75, cex.key.xlab = 0.75,
    cex.key.ylab = 0.75,
    density.info = c("histogram", "density", "none"),
    denscol = tracecol, densadj = 0.25, main = "Heatmap",
    sub = "", xlab = "", ylab = "", cex.main = 2,
    cex.sub = 1.5, font.main = 2, font.sub = 3,
    adj.main = 0.5, mgp.main = c(1.5, 0.5, 0),
    mar.main = 3, mar.sub = 3, if.plot = TRUE,
    plot.row.partition = FALSE, plot.col.partition = FALSE,
    cex.partition = 1.25, color.partition.box = "gray45",
    color.partition.border = "#FFFFFF",
    plot.row.individuals = FALSE,
    plot.col.individuals = FALSE,
    plot.row.clusters = FALSE, plot.col.clusters = FALSE,
    plot.row.clustering = FALSE,
    plot.col.clustering = FALSE,
    plot.row.individuals.list = FALSE,
    plot.col.individuals.list = FALSE,
    plot.row.clusters.list = FALSE,
    plot.col.clusters.list = FALSE,
    plot.row.clustering.list = FALSE,
    plot.col.clustering.list = FALSE, row.data = FALSE,
    col.data = FALSE, if.plot.info = FALSE, text.box,
    cex.text = 1, ...)
}
\arguments{
  \item{x}{data matrix or data frame, or dissimilarity
  matrix or `dist' object determined by the value of the
  'diss' argument. ##diss logical flag: if TRUE (default
  for \code{dist} or \code{dissimilarity} objects), then
  \code{x} is assumed to be a dissimilarity matrix. If
  FALSE,then \code{x} is treated as a matrix of
  observations by variables.}

  \item{diss}{logical, whether the \code{x} is a
  dissimilarity matrix}

  \item{Rowv}{one of the following: TRUE, a `dend' object,
  a vector or NULL/FALSE; determines if and how the
  \emph{row} dendrogram should be reordered.}

  \item{Colv}{one of the following: "Rowv", TRUE, a `dend'
  object, a vector or NULL/FALSE; determines if and how the
  \emph{column} dendrogram should be reordered.}

  \item{dendrogram}{character string indicating whether to
  draw 'none', 'row', 'column' or 'both' dendrograms.
  Defaults to 'both'.}

  \item{dist.row}{a \code{dist} object for \emph{row}
  observations.}

  \item{dist.col}{a \code{dist} object for \emph{column}
  observations.}

  \item{dist.FUN}{function used to compute the distance
  (dissimilarity) between both rows and columns.  Defaults
  to \code{gdist}.}

  \item{dist.FUN.MoreArgs}{a list of other arguments to be
  passed to \code{gdist}}

  \item{hclust.row}{a \code{hclust} object (as produced by
  \code{hclust}) for \emph{row} observations.}

  \item{hclust.col}{a \code{hclust} object (as produced by
  \code{hclust}) for \emph{column} observations.}

  \item{hclust.FUN}{function used to compute the
  hierarchical clustering when "Rowv" or "Colv" are not
  dendrograms. Defaults to \code{hclust}.}

  \item{hclust.FUN.MoreArgs}{a list of other arguments to
  be passed to \code{hclust}.\cr Defaults to
  \code{list(method="ward")}}

  \item{scale}{character indicating if the values should be
  centered and scaled in either the row direction or the
  column direction, or none.  The default is
  \code{"none"}.}

  \item{na.rm}{logical, whether NA values will be removed
  when scaling.}

  \item{cluster.by.row}{logical, whether to cluster
  \emph{row} observations and reorder the input
  accordingly.}

  \item{cluster.by.col}{logical, whether to cluster
  \emph{column} observations and reorder the input
  accordingly.}

  \item{kr}{numeric, number of clusters in rows; suppressed
  when \code{row.cluster} is specified. DEFAULT: NULL.}

  \item{kc}{numeric, number of clusters in columns;
  suppressed when \code{col.cluster} is specified. DEFAULT:
  NULL.}

  \item{row.clusters}{a numerical vector, indicating the
  cluster labels of \emph{row} observations.}

  \item{col.clusters}{a numerical vector, indicating the
  cluster labels of \emph{column} observations.}

  \item{revR}{logical indicating if the row order should be
  'rev'ersed for plotting.}

  \item{revC}{logical indicating if the column order should
  be 'rev'ersed for plotting, such that e.g., for the
  symmetric case, the symmetry axis is as usual.}

  \item{add.expr}{expression that will be evaluated after
  the call to \code{image}. Can be used to add components
  to the plot.}

  \item{breaks}{numeric, either a numeric vector indicating
  the splitting points for binning \code{x} into colors, or
  a integer number of break points to be used, in which
  case the break points will be spaced equally between
  \code{range(x)}. DEFAULT: 16 when not specified.}

  \item{x.center}{numeric, a value of \code{x} for
  centering colors to}

  \item{color.FUN}{function or function name in characters,
  for colors in the heatmap}

  \item{sepList}{a \code{list} of length 2 giving the row
  and column lines of separation.}

  \item{sep.color}{color for lines of separation.}

  \item{sep.lty}{line type for lines of separation.}

  \item{sep.lwd}{line width for lines of separation.}

  \item{cellnote}{(optional) matrix of character strings
  which will be placed within each color cell, e.g. cell
  labels or p-value symbols.}

  \item{cex.note}{relative font size of \code{cellnote}.}

  \item{notecol}{color of \code{cellnote}.}

  \item{na.color}{Color to use for missing value
  (\code{NA}). Defaults to the plot background color.}

  \item{trace}{character string indicating whether a solid
  "trace" line should be drawn across \code{"row"}s or down
  \code{"column"}s, \code{"both"} or \code{"none"}. The
  distance of the line from the center of each color-cell
  is proportional to the size of the measurement. Defaults
  to \code{"none"}.}

  \item{tracecol}{character string giving the color for
  "trace" line. Defaults to "cyan";}

  \item{hline}{Vector of values within cells where a
  horizontal dotted line should be drawn. only plotted if
  'trace' is 'row' or 'both'. Default to the median of the
  breaks.}

  \item{vline}{Vector of values within cells where a
  vertical dotted line should be drawn; only drawn if
  'trace' 'column' or 'both'. \code{vline} default to the
  median of the breaks.}

  \item{linecol}{the color of \code{hline} and
  \code{vline}. Defaults to the value of 'tracecol'.}

  \item{labRow}{character vectors with row labels to use;
  defaults to \code{rownames(x)}.}

  \item{labCol}{character vectors with column labels to
  use; defaults to \code{colnames(x)}.}

  \item{margin.for.labRow}{a numerical value gives the
  margin to plot \code{labRow}.}

  \item{margin.for.labCol}{a numerical value gives the
  margin to plot \code{labCol}.}

  \item{ColIndividualColors}{(optional) character vector of
  length \code{ncol(x)} containing the color names for a
  horizontal side bar that may be used to annotate the
  columns of \code{x}.}

  \item{RowIndividualColors}{(optional) character vector of
  length \code{nrow(x)} containing the color names for a
  vertical side bar that may be used to annotate the rows
  of \code{x}.}

  \item{cexRow}{positive numbers, used as 'cex.axis' in for
  column axis labeling. The default currently only uses
  number of columns.}

  \item{cexCol}{positive numbers, used as 'cex.axis' in for
  the row axis labeling. The default currently only uses
  number of rows.}

  \item{labRow.by.group}{logical, whether group unique
  labels for rows.}

  \item{labCol.by.group}{logical, whether group unique
  labels for columns.}

  \item{key}{logical indicating whether a color-key should
  be shown.}

  \item{keysize}{numeric value indicating the relative size
  of the key}

  \item{mapsize}{numeric value indicating the relative size
  of the heatmap.}

  \item{mapratio}{the width-to-height ratio of the
  heatmap.}

  \item{sidesize}{numeric value indicating the relative
  size of the sidebars.}

  \item{cex.key.main}{a numerical value giving the amount
  by which \code{main}-title of color-key should be
  magnified relative to the default.}

  \item{cex.key.xlab}{a numerical value giving the amount
  by which \code{xlab} of color-key should be magnified
  relative to the default.}

  \item{cex.key.ylab}{a numerical value giving the amount
  by which \code{ylab} of color-key should be magnified
  relative to the default.}

  \item{density.info}{character string indicating whether
  to superimpose a 'histogram', a 'density' plot, or no
  plot ('none') on the color-key.}

  \item{denscol}{character string giving the color for the
  density display specified by 'density.info', defaults to
  the same value as 'tracecol'.}

  \item{densadj}{Numeric scaling value for tuning the
  kernel width when a density plot is drawn on the color
  key.  (See the 'adjust' parameter for the 'density'
  function for details.)  Defaults to 0.25.}

  \item{main}{an overall title for the plot. See
  \code{help("title", package="graphics")}.}

  \item{sub}{a subtitle for the plot, describing the
  distance and/or alignment gap (the "shift").}

  \item{xlab}{a title for the x axis. See
  \code{help("title", package="graphics")}.}

  \item{ylab}{a title for the y axis. See
  \code{help("title", package="graphics")}.}

  \item{cex.main}{a numerical value giving the amount by
  which \code{main}-title should be magnified relative to
  the default.}

  \item{cex.sub}{a numerical value giving the amount by
  which \code{sub}-title should be magnified relative to
  the default.}

  \item{font.main}{An integer which specifies which font to
  use for \code{main}-title.}

  \item{font.sub}{An integer which specifies which font to
  use for \code{sub}-title.}

  \item{adj.main}{The value of 'adj' determines the way in
  which \code{main}-title strings are justified.}

  \item{mgp.main}{the margin line (in 'mex' units) for the
  \code{main}-title.}

  \item{mar.main}{a numerical vector of the form
  \code{c(bottom, left, top, right)} which gives the number
  of lines of margin to be specified on the four sides of
  the \code{main}-title.}

  \item{mar.sub}{a numerical vector of the form
  \code{c(bottom, left, top, right)} which gives the number
  of lines of margin to be specified on the four sides of
  the \code{sub}-title.}

  \item{if.plot}{logical, whether to plot. Reordered matrix
  is returned without graphical output if FALSE.}

  \item{plot.row.partition}{logical, whether to plot
  \emph{row} partition.}

  \item{plot.col.partition}{logical, whether to plot
  \emph{column} partition.}

  \item{cex.partition}{a numerical value giving the amount
  by which \code{partition} should be magnified relative to
  the default.}

  \item{color.partition.box}{color for the \code{partition}
  box.}

  \item{color.partition.border}{color for the
  \code{partition} border.}

  \item{plot.row.individuals}{logical, whether to make a
  plot of \emph{row} observations.}

  \item{plot.col.individuals}{logical, whether to make a
  plot of \emph{column} observations.}

  \item{plot.row.clusters}{logical, whether to make a
  summary plot of \emph{row} clusters.}

  \item{plot.col.clusters}{logical, whether to make a
  summary plot of \emph{column} clusters.}

  \item{plot.row.clustering}{logical, whether to make a
  summary plot of overall \emph{row} clustering.}

  \item{plot.col.clustering}{logical, whether to make a
  summary plot of overall \emph{column} clustering.}

  \item{plot.row.individuals.list}{a list of expressions
  that is used to \code{plot.row.individuals}}

  \item{plot.col.individuals.list}{a list of expressions
  that is used to \code{plot.col.individuals}}

  \item{plot.row.clusters.list}{a list of expressions that
  is used to \code{ plot.row.clusters}}

  \item{plot.col.clusters.list}{a list of expressions that
  is used to \code{plot.col.clusters}}

  \item{plot.row.clustering.list}{a list of expressions
  that is used to \code{plot.row.clustering}}

  \item{plot.col.clustering.list}{a list of expressions
  that is used to \code{plot.col.clustering}}

  \item{row.data}{(optional) data used to
  \code{plot.row.individuals}, \code{ plot.row.clusters} or
  \code{plot.row.clustering}}

  \item{col.data}{(optional) data used to
  \code{plot.col.individuals}, \code{ plot.col.clusters} or
  \code{plot.col.clustering}}

  \item{if.plot.info}{logical, whether to plot
  \code{text.box}.}

  \item{text.box}{character plotted when
  \code{if.plot.info} is TRUE.}

  \item{cex.text}{a numerical value giving the amount by
  which \code{text.box} should be magnified relative to the
  default.}

  \item{...}{arguments to be passed to method
  \code{heatmap.3}.\cr e \code{help("image",
  package="graphics")}.}
}
\value{
  A reordered matrix according to \emph{row} or/and
  \emph{col} dendrogram(s) and indices that used for
  reordering.
}
\description{
  Enhanced heatmap representation with dendrograms and
  partition given the \emph{elbow criterion} or a desired
  number of clusters.\cr 1) a dendrogram added to the left
  side and to the top, according to cluster analysis;\cr 2)
  partitions in highlighted rectangles, according to the
  "elbow" rule or a desired number of clusters.\cr
}
\details{
  Enhanced heatmap representation with partition and
  summary statistics (optional). This is an enhanced
  version of `heatmap.2' function in the Package
  \pkg{gplots}. The enhancement includes: 1) Improved
  performance with optional input of precomputed
  \emph{\code{dist}} object and \emph{\code{hclust}}
  object. 2) Highlight of specific cells using rectangles.
  For instance, the cells of clusters of interests.
  (Examples should be included in future.) 3) Add-on plots
  in addition to the heatmap, such as cluster-wise summary
  plots and overall clustering summary plots, to the right
  of or under the heatmap.
}
\examples{
## ------------------------------------------------------------------------
## Example1: mtcars
## ------------------------------------------------------------------------
## load library
require("GMD")

## load data
data(mtcars)

## heatmap on raw data
x  <- as.matrix(mtcars)

dev.new(width=10,height=8)
heatmap.3(x)                               # default, with reordering and dendrogram
heatmap.3(x, Rowv=FALSE, Colv=FALSE)       # no reordering and no dendrogram
heatmap.3(x, dendrogram="none")            # reordering without dendrogram
heatmap.3(x, dendrogram="row")        # row dendrogram with row (and col) reordering
heatmap.3(x, dendrogram="row", Colv=FALSE) # row dendrogram with only row reordering
heatmap.3(x, dendrogram="col")             # col dendrogram
heatmap.3(x, dendrogram="col", Rowv=FALSE) # col dendrogram with only col reordering
heatmapOut <-
  heatmap.3(x, scale="column")             # sacled by column
names(heatmapOut)                          # view the list that is returned
heatmap.3(x, scale="column", x.center=0)   # colors centered around 0
heatmap.3(x, scale="column",trace="column")  # trun "trace" on

## coloring cars (row observations) by brand
brands <- sapply(rownames(x), function(e) strsplit(e," ")[[1]][1])
names(brands) <- c()
brands.index <- as.numeric(as.factor(brands))
RowIndividualColors <- rainbow(max(brands.index))[brands.index]
heatmap.3(x, scale="column", RowIndividualColors=RowIndividualColors)

## coloring attributes (column features) randomly (just for a test :)
heatmap.3(x, scale="column", ColIndividualColors=rainbow(ncol(x)))

## add a single plot for all row individuals
dev.new(width=12,height=8)
expr1 <- list(quote(plot(row.data[rowInd,"hp"],1:nrow(row.data),
xlab="hp",ylab="",yaxt="n",main="Gross horsepower")),
quote(axis(2,1:nrow(row.data),rownames(row.data)[rowInd],las=2)))
heatmap.3(x, scale="column", plot.row.individuals=TRUE, row.data=x,
          plot.row.individuals.list=list(expr1))


## ------------------------------------------------------------------------
## Example2: ruspini
## ------------------------------------------------------------------------
## load library
require("GMD")
require(cluster)

## load data
data(ruspini)

## heatmap on a `dist' object
x <- gdist(ruspini)
main <- "Heatmap of Ruspini data"
dev.new(width=10,height=10)
heatmap.3(x, main=main, mapratio=1) # with a title and a map in square!
heatmap.3(x, main=main, revC=TRUE)  # reverse column for a symmetric look
heatmap.3(x, main=main, kr=2, kc=2) # partition by predefined number of clusters

## show partition by elbow
css.multi.obj <- css.hclust(x,hclust(x))
elbow.obj <- elbow.batch(css.multi.obj,ev.thres=0.90,inc.thres=0.05)
heatmap.3(x, main=main, revC=TRUE, kr=elbow.obj$k, kc=elbow.obj$k)

## show elbow info as subtitle
heatmap.3(x, main=main, sub=sub("\\n"," ",attr(elbow.obj,"description")),
cex.sub=1.25,revC=TRUE,kr=elbow.obj$k, kc=elbow.obj$k)
}

