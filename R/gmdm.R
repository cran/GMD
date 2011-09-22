
##' Generalized Minimum Distance Matrix
##'
##' Generalized Minimum Distance Matrix
##' @title Generalized Minimum Distance Matrix
##' @param x a list of numeric vectors
##' @param labels a character vector of the same length of x, giving the names of the numeric vectors.
##' @param pseudocount a numeric value to be allocated for each position to reduce bias;
##' by default \code{pseudocount = 0}.
##' @param sliding logical, indicating whether sliding is allowed or not for an optimal solution; 
##' by default \code{sliding = TRUE}.
##' @return \code{gmdm} returns an object of class \code{gmdm}, a list with components
##'
##' labels: a string vector, giving the names of distributions
##'
##' data.ori: a list of numeric vectors, giving the original input
##'
##' data: a list of numeric vectors, giving the normalized version of the original input
##' 
##' dm: a numeric numeric, the pairwise distance matrix of \emph{GM-Distances}
##' 
##' gap.pair: a numeric matrix, giving the gap pair of each alignment per row:
##' i.e. relative shifts between distributions of the optimal hit
##' 
##' sliding: logical, indicating whether sliding is performed
##' 
##' pseudocount: a numeric value that is allocated at each position in addition to original values
##' @references See \code{citation("GMD")}
##' @seealso \code{\link{plot.gmdm}}, \code{\link{gmd}}
##' @keywords classes
##' @examples
##' require(GMD)
##' data(cage)
##' x <- gmdm(cage)
##' print(x$labels)
##' print(x$dm)
##' 
##' \dontrun{data(cagel)
##' x <- gmdm(cagel)
##' head(x$labels)
##' head(x$dm)}
gmdm <- 
  function(x,
           labels=names(x),
           pseudocount=0,
           sliding=TRUE
           )
{
  ret <- list()
  N <- length(x)
  dm <- matrix(0, ncol=N, nrow=N)
  colnames(dm) <- rownames(dm) <- labels
  gap.pair <- matrix(0, ncol=2, nrow=N*(N-1)/2)
  n <- 0
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      n <- n+1
      ij.x <- gmd(x[[i]],x[[j]],pseudocount=pseudocount,sliding=sliding)
      dm[i,j] <- dm[j,i] <- ij.x$d
      gap.pair[n,] <- ij.x$gap.pair[1,]
      if (nrow(ij.x$gap.pair)>1){
        warning(sprintf("There are multiple (equally good) optimal hits/alignment
between %s and %s;
the first one will be saved.",labels[i],labels[j]))
      }
    }
  }
  ret$labels <- labels
  ret$data.ori <- x
  ret$data <- lapply(x,function(e) e/sum(e))
  ret$dm <- dm
  ret$gap.pair <- gap.pair
  ret$sliding <- sliding
  ret$pseudocount <- pseudocount
  class(ret) <- c("gmdm","list")
  return(ret)
}







##' Plot Function for Class \code{gmdm}
##'
##' Plot Function for Class \code{gmdm}
##' @title Plot Function for Class gmdm
##' @param x an object of class \code{gmdm}.
##' @param labels a string vector of the same length as \code{x$data},
##' giving the names of the numeric vectors in \code{x$data}.
##' @param colors the colors of the discrete distributions; the default is \emph{"Dark2" colors in ColorBrewer palettes} if not specified.
##' @param plot.type the plot type. See \code{type} in \code{plot} for possible values; the default \code{plot.type = "h"}, giving \sQuote{\bold{h}istogram} like vertical lines.
##' @param main an overall title for the plot. See \code{help("title", package="graphics")}; the default title is used if not specified.
##' @param ylab a title for the y axis. See \code{help("title", package="graphics")}.
##' @param xlab a title for the x axis. See \code{help("title", package="graphics")}.
##' @param label.length.max numeric, giving the maximum string width allowed in diagonal labels.
##' @param label.line.max numeric, giving the maximum number of lines allowed in diagonal labels.
##' @param cex.text a numerical value giving the amount by which plot text should be magnified relative to the default.
##' @param cex.tickmark a numerical value giving the amount by which tickmarks should be magnified relative to the default.
##' @param if.plot.new logical, indicating whether to start a new plot device.
##' @param x.jitter numeric, indicating how \emph{jitter} should be added to distinguish subplots; by default \code{x.jitter=1/1000} indicating the jitter is adjusted to 1/1000 of the \emph{x-axis range}.
##' @param ... arguments to be passed to methods, see \code{gmd}.
##' @references See \code{help(GMD)}
##' @seealso \code{\link{gmdm}}, \code{\link{gmd}}
##' @keywords methods hplot
##' @examples
##' require(GMD)
##' data(cage)
##' plot(gmdm(cage))
plot.gmdm <- 
  function (x,
            ##
            labels=x$labels,
            colors,
            ##
            plot.type="h",
            ##
            main,
            ylab="Fraction",
            xlab="Position",
            ## maximum size of diagonal labels
            label.length.max=8,
            label.line.max=3,
            ## style
            cex.text=2,
            cex.tickmark=0.75,
            ##
            if.plot.new=TRUE,
            x.jitter=1/1000,
            ...
            ) 
{
  if (!.is.gmdm(x)){
    stop("`x' should be an object of class `gmdm'.")
  }
  
  ## params ##
  if(missing(colors)){
    colors <- .brewer.pal.Dark2 
  }
  
  N <- nrow(x$dm)
  
  ## labels for plot ##
  nchar.max <- max(nchar(labels))
  label.size.max <- label.length.max*label.line.max
  if(nchar.max>label.size.max){
    warning(sprintf("The number of characters in each `label` should not exceed %s.
The first %s characters are kept.",label.size.max)
            )
  }
  labels.plot <- sapply(labels,substr,start=1,stop=label.size.max)
  labels.plot <- sapply(labels.plot,.wordwrap,len=label.length.max)

  ## colors ##
  if(length(colors) < N){
    warning("The length of `colors' is shorter than that of `x' and the `colors' are recycled.")
    ##tmp <- cbind(colors,1:N)
    colors <- rep(colors,ceiling(N/length(colors)))
  }

  ## Main title ##
  if(missing(main)){
    main <-
      sprintf("Optimal alignments among distributions (%s sliding)",
              ifelse(x$sliding,"with","without")
              )
  }
  
  ## plot ## 
  ylim.max <- max(unlist(lapply(x$data,max)))
  xlim.max <- max(unlist(lapply(x$data,length)))
  ## print(sprintf("ylim.max=%s",ylim.max))
  ## print(sprintf("xlim.max=%s",xlim.max))

  ## 
  if (if.plot.new) {
    dev.new(width=8,height=8)
  }
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  par(mfrow=c(N,N),
      mgp=c(0,0.2,0),
      mar=c(1,1,1,1),
      oma=c(2,2,4,1))

  n <- 0
  for (i in 1:(N-1)){
    par(mfg=c(i,i))
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
    text(1,1,labels.plot[i],cex=cex.text/N*3.6*1.2,col=colors[i])
    for (j in (i+1):N) {
      n <- n+1
      par(mfg=c(i,j))
      ij.x <- list()
      ij.x$labels <- x$labels[c(i,j)]
      ij.x$v1.ori <- x$data.ori[[i]]
      ij.x$v2.ori <- x$data.ori[[j]]
      ij.x$v1 <- x$data[[i]]
      ij.x$v2 <- x$data[[j]]
      ij.x$distance <- x$dm[i,j]
      ij.x$sliding <- x$sliding
      ij.x$pseudocount <- x$pseudocount
      ij.x$gap.pair <- matrix(x$gap.pair[n,],nrow=1)
      ij.x$n.hit <- 1
      class(ij.x) <- c("gmd","list")
      
      plot.gmd(x=ij.x,
           labels=ij.x$labels,
           colors=colors[c(i,j)],
           plot.method="overlay",
           main="",
           if.plot.new=FALSE,
           if.text.gmd=FALSE,
           if.text.gap=FALSE,
           if.plot.gap=FALSE,
           if.plot.legned=FALSE,
           xlab="",
           ylab="",
           ylim=c(0,ylim.max),
           xlim=c(1,xlim.max),
           cex.tickmark=cex.tickmark,
           ...
           )
      ##
      par(mfg=c(j,i))
      ij.text <-
        sprintf("GMD=%.3f\nGap=%s",ij.x$distance,sprintf("(%s,%s)",ij.x$gap.pair[1],ij.x$gap.pair[2]))
      plot(c(0,1),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n")
      text(0.5,0.5,ij.text,cex=cex.text/N*3.6*1.05)
    }
  }
  par(mfg=c(N,N))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  text(1,1,labels.plot[N],cex=cex.text/N*3.6*1.2,col=colors[N])
  mtext(ylab,line=0,side=2,outer=TRUE,cex=1.5)
  mtext(xlab,line=0,side=1,outer=TRUE,cex=1.5)
  mtext(text=main,line=1.0,side=3,outer=TRUE,cex=2,adj=0.5)

  invisible()
  
}
