
## ************************************************************************
##
## Compute GM-Distance between discrete distributions, using method from:
##
## Systematic Clustering of Transcription Start Site Landscapes
## Author(s): Zhao et al, 2011
## Source: PLoS ONE 6(8): e23409. doi:10.1371/journal.pone.0023409 
## 
## URL: http://dx.plos.org/10.1371/journal.pone.0023409
##      http://www.plosone.org/article/info:doi/10.1371/journal.pone.0023409
##
##
## (c) Xiaobei Zhao 2011 <xiaobei@binf.ku.dk>
## 
##require(tools)
## ##RColorBrewer
## ************************************************************************

##' Generalized Minimum Distance
##'
##' Generalized Minimum Distance
##' @title Generalized Minimum Distance (GMD)
##' @param v1 a numeric vector, giving positional counts as a discrete distribution.
##' @param v2 a numeric vector, giving positional counts as a discrete distribution.
##' @param labels a string vector of length 2, giving the names of v1 and v2 respectively.
##' @param pseudocount a numeric value to be allocated for each position to reduce bias;
##' by default \code{pseudocount = 0}.
##' @param sliding logical, indicating whether sliding is allowed or not for an optimal solution; 
##' by default \code{sliding = TRUE}.
##' @return \code{gmd} returns an object of class \code{gmd}, a list with components
##'
##' labels: a string vector, giving the names of distributions
##' 
##' v1.ori: a numeric vector, the first input distribution
##' 
##' v2.ori: a numeric vector, the second input distribution
##' 
##' v1: a numeric vector, the normalized version of the first input distribution
##' 
##' v2: a numeric vector, the normalized version of the second input distribution
##' 
##' distance: numeric, the \emph{GM-Distance} (\emph{GMD})
##' 
##' sliding: logical, indicating whether sliding is performed
##' 
##' pseudocount: a numeric value that is allocated at each position in addition to original values
##' 
##' gap.pair: a numeric matrix, giving one gap pair per row:
##' i.e. relative shifts between distributions of one optimal hit
##' 
##' n.hit: numeric, the number of (equally good) optimal hits
##' @references See \code{citation("GMD")}
##' @seealso \code{\link{print.gmd}}, \code{\link{summary.gmd}}, \code{\link{plot.gmd}}, \code{\link{gmdm}} 
##' @keywords classes
##' @examples
##' require(GMD)
##' gmd(c(4,1,1,0,0,0,3,1),c(2,1,1,0,0,0,3,3),sliding=FALSE)
##' x <- gmd(c(4,1,1,0,0,0,3,1), c(1,1,2,1,1,0,0,0,3,3,5,5),
##' pseudocount=1, sliding=TRUE)
##' print(x)
##' print(x, "full")
gmd <-
  function(v1, v2, labels=c("v1","v2"), pseudocount=0, sliding=TRUE)
{
  ## 
  ret <- list()
  ret$labels <- labels
  ret$v1.ori <- v1
  ret$v2.ori <- v2
  
  ##
  l1 <- length(v1)
  l2 <- length(v2)

  ##
  if (sliding) {
    res <- .gmd(v1, v2, pseudocount)
    sliding <- TRUE
    gap <- which(res$position==1) - min(l1,l2) - 1
    n.hit <- sum(res$position==1)
  } else {
    if (l1 != l2){
      stop("The lengths of two input vectors should be equal.
Otherwise, please allow `sliding` (sliding=TRUE).\n"
           )
    } else{
      d <- .gmd0(v1, v2, pseudocount)
      p <- 0
      res <- list(distance=d, position=c(1))
      sliding <- FALSE
      gap <- 0
      n.hit <- 1 
    }
  }

  if (l1<=l2){
    gap.pair <- cbind(gap,0)
  } else {
    gap.pair <- cbind(0,gap)
  }
  colnames(gap.pair) <- c()
  gap.pair <- t(apply(gap.pair,1,function(x) if(any(x<0)) x-min(x) else x))

  ret$v1 <- v1/sum(v1)
  ret$v2 <- v2/sum(v2)
  ret$distance <- res$d
  ret$sliding <- sliding
  ret$pseudocount <- pseudocount
  ret$gap.pair <- gap.pair
  ret$n.hit <- n.hit
  class(ret) <- c("gmd","list")
  return(ret)
}



##' Print Function for Class \code{gmd}
##'
##' Print Function for Class \code{gmd}
##' @title Print Function for Class gmd
##' @param x an object of class \code{gmd}.
##' @param print.mode a string, indicating whether to print in \emph{full} mode (\emph{default}). 
##' @param digits integer, indicating the number of decimal places to be printed.
##' @param ... arguments to be passed to methods, see \code{print}.
##' @seealso \code{\link{gmd}}
##' @references See \code{help(GMD)}
print.gmd <- function(x, print.mode=c("brief","full"), digits=3, ...)
{
  if (!.is.gmd(x)){
    stop("`x' should be an object of class `gmd'.")
  }
  print.mode <- match.arg(print.mode)
  s.gap <-
    paste(sprintf("\t%s\t%s\n",x$labels[1],x$labels[2]),
          paste(sapply(1:nrow(x$gap.pair),
                       function(i) sprintf("Hit%s\t%s\t%s",i, x$gap.pair[i,1],x$gap.pair[i,2])),
                sep="",
                collapse="\n"),
          "\n",
          sep=""
          )
  cat("\n")
  if(print.mode=="full"){
    cat(sprintf("Distribution of %s:\n%s\n",x$labels[1],paste(round(x$v1.ori,digits),sep="",collapse=" ")))
    cat(sprintf("(After normalization)\n%s\n\n",paste(round(x$v1,digits),sep="",collapse=" ")))
    cat(sprintf("Distribution of %s:\n%s\n",x$labels[2],paste(round(x$v2.ori,digits),sep="",collapse=" ")))
    cat(sprintf("(After normalization)\n%s\n\n",paste(round(x$v2,digits),sep="",collapse=" ")))
  }
  cat(sprintf(sprintf("GM-Distance: %%.%sf\n\n",digits),x$distance))
  cat(sprintf("Sliding: %s\n\n",x$sliding))
  cat(sprintf("Number of hits: %s\n\n",x$n.hit))
  cat(sprintf("Gap:\n%s\n\n",s.gap))
  cat("\n")
}


 
##' Summary Function for Class \code{gmd}
##'
##' Summary Function for Class \code{gmd}
##' @title Summary Function for Class gmd
##' @param object an object of class \code{gmd}.
##' @param ... arguments to be passed to methods, see \code{summary}.
##' @seealso \code{\link{gmd}}
##' @references See \code{help(GMD)}
summary.gmd <- function(object, ...){
  if (!.is.gmd(object)){
    stop("`object' should be an object of class `gmd'.")
  }
  class(object) <- "list"
  summary(object, ...)
}



##' Plot Function for Class \code{gmd}
##'
##' Plot Function for Class \code{gmd}
##' @title Plot Function for Class gmd
##' @param x an object of class \code{gmd}.
##' @param labels a string vector of the same length of \code{x$labels},
##' giving the names of the numeric vectors in \code{x}.
##' @param colors the colors of the discrete distributions; by default they are in \code{"red"} and \code{"blue"}.
##' @param plot.method the plot method. This can be specified as a string:
##' \code{"separate"}: means separated subplots [\emph{default}];
##' \code{"overlay"}: means overlaid subplots.
##' @param plot.type the plot type. See \code{type} in \code{plot} for possible values; the default \code{plot.type = "h"}, giving \sQuote{\bold{h}istogram} like vertical lines.
##' @param main an overall title for the plot. See \code{help("title", package="graphics")}.
##' @param ylab a title for the y axis. See \code{help("title", package="graphics")}.
##' @param xlab a title for the x axis. See \code{help("title", package="graphics")}.
##' @param ylim range of y values, as in \code{help("plot", package="graphics")}.
##' @param xlim range of x values, as in \code{help("plot", package="graphics")}.
##' @param font.type the name of a font type for drawing text. See \code{font} in \code{par}; the default \code{font.type = 1}, corresponding to plain text.
##' @param font.family the name of a font family for drawing text. See \code{family} in \code{par}; the default \code{font.family = "sans"}, corresponding to san serif typeface.
##' @param cex.lab a numerical value giving the amount by which \code{xlab} and \code{ylab} should be magnified relative to the default.
##' @param cex.tickmark a numerical value giving the amount by which tickmarks should be magnified relative to the default.
##' @param cex.legend a numerical value giving the amount by which legends should be magnified relative to the default.
##' @param lwd.line the line width, a \emph{positive} number, defaulting to \code{1}.
##' @param if.plot.new logical, indicating whether to start a new plot device.
##' @param if.text.gmd logical, indicating whether \emph{GM-Distance} is reported in the subtitle.
##' @param if.text.gap logical, indicating whether \emph{gap} is reported in the subtitle.
##' @param if.plot.gap logical, indicating whether \emph{gap} is plotted.
##' @param if.plot.legned logical, indicating whether \emph{legend} is plotted.
##' @param x.jitter numeric, indicating how \emph{jitter} should be added to distinguish subplots; by default \code{x.jitter=ifelse(plot.method=="overlay",1/1000,0)} giving how jitter should be adjusted according to the \emph{x-axis range}.
##' @param ... arguments to be passed to methods, such as graphical parameters (see \code{par}).
##' @references See \code{help(GMD)}
##' @seealso \code{\link{gmd}}
##' @keywords methods hplot
##' @examples
##' require(GMD)
##' data(cage)
##' \dontrun{plot(gmd(cage[[1]],cage[[2]],labels=names(cage)[c(1,2)],
##' pseudocount=1, sliding=TRUE))}
##' plot(gmd(cage[[1]],cage[[3]],labels=names(cage)[c(1,3)],
##' pseudocount=1, sliding=TRUE))
##' plot(gmd(cage[[1]],cage[[3]],labels=names(cage)[c(1,3)],
##' pseudocount=1, sliding=TRUE), plot.method="overlay")
plot.gmd <-
  function(x,
           ##
           labels=x$labels,
           colors=c("red","blue"),
           ##
           plot.method=c("separate","overlay"),
           plot.type="h",
           ##
           main,
           ylab="Fraction",
           xlab="Position",
           ylim, 
           xlim,
           ## style
           font.type=1,
           font.family=c("sans","serif","mono"),
           cex.lab=1.2,
           cex.tickmark=1,
           cex.legend=1.5,
           lwd.line=1,
           ##
           if.plot.new=TRUE,
           if.text.gmd=TRUE,
           if.text.gap=FALSE,
           if.plot.gap=TRUE,
           if.plot.legned=TRUE,
           ##
           x.jitter=ifelse(plot.method=="overlay",1/1000,0),
           ...)
{

  ## --
  ## Top: v1
  ## Middle: alignment
  ## Bottom: v2
  ## --
  if (!.is.gmd(x)){
    stop("`x' should be an object of class `gmd'.")
  }
  
  ## params ##
  font.family <- match.arg(font.family)
  plot.method <- match.arg(plot.method)

  ## check multiple hits ##
  if (nrow(x$gap.pair)>1){
    warning("There are multiple optimal hits/alignment; the first one will be plotted.")
  }
  gap.pair <- x$gap.pair[1,] ## pick the first if multiple cases
  s.gap <- sprintf("gap=%s",sprintf("(%s,%s)",gap.pair[1],gap.pair[2]))

  ## title and subtitle ##
  if(missing(main)){
    main <-
      sprintf("Optimal alignment between distributions (%s sliding)",
              ifelse(x$sliding,"with","without")
              )
  }
  
  s.gmd <- sprintf("GMD=%.3f",x$distance)
  
  if (if.text.gmd & if.text.gap){
    s.sub <- sprintf("%s (%s)",s.gmd,s.gap)
  } else if (if.text.gmd){
    s.sub <- s.gmd
  } else if (if.text.gap){
    s.sub <- s.gap
  } else {
    s.sub <- ""
  }

  
  
  ## axes ##
  if (missing(xlim)){
    xlim <- c(min(1,length(x$v2)+x$gap),max(length(x$v2)+x$gap,length(x$v1)))
  }
  if (missing(ylim)){
    ylim <- c(0, max(x$v2, x$v1))
  }
  x.breaks <- unique(c(1,pretty(xlim, n=4, high.u.bias=10)))
  x.breaks <- x.breaks[x.breaks>=1]
  y.breaks <- pretty(ylim,n=4,high.u.bias=10)

  
  ## alignment ##
  s.align <- c(rep("-",max(gap.pair)),rep(" ",xlim[2]-xlim[1]+1-max(gap.pair[1])))
  
  
  ## plot
  if (if.plot.new) {
    dev.new(width=10, height=8)
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(font.main=font.type, font.lab=font.type, font.axis=font.type)
    par(family=font.family)
  }
  
  if (plot.method=="separate"){
    ## separate subplot
    if (if.plot.new) {
      par(mfrow=c(2,1),mar=c(1.5,4,1.8,0),omi=c(0.5,0.25,0.75,0.5),mgp=c(2.5,0.25,0))
    }
    ##1
    plot(x$v1,
         type=plot.type,
         ylim=ylim,
         xlim=xlim-gap.pair[1],
         ylab="",
         yaxt="n",
         xlab="",
         xaxt="n",
         col=colors[1],
         lwd=lwd.line
         )
   
    axis(side=2,at=y.breaks,labels=y.breaks,tck=-0.015,cex.axis=cex.tickmark,las=1)  # y-axis
    axis(side=1,at=x.breaks+gap.pair[2],labels=x.breaks,tck=-0.015,cex.axis=cex.tickmark) # x-axis
    axis(1, at=x.breaks,labels=rep("",length(x.breaks)),tck=-0.015) # alignment
    mtext(line=2,side=2,text=ylab,cex=cex.lab) # y-lab
    
    mtext(line=2,side=3,text=main,outer=FALSE,cex=2,adj=0) # main title
    mtext(line=0.55,side=3,text=s.sub,outer=FALSE,cex=1.5,adj=0) # subtitle
    
    if(if.plot.legned)
      legend("topright",legend=labels[1],text.col=colors[1],cex=cex.legend,inset=0,box.col=FALSE) # legend
    if (if.plot.gap){
      mtext(line=1.5,side=1,at=seq(xlim[1],xlim[2])-gap.pair[1],s.align,cex=3) # gap/alignment
    }
    
    ##2
    plot(x$v2,
         type=plot.type,
         ylim=ylim,
         xlim=xlim-gap.pair[2],
         ylab="",
         yaxt="n",
         xlab="",
         xaxt="n",
         col=colors[2],
         lwd=lwd.line
         )
    
    axis(side=2,at=y.breaks,labels=y.breaks,tck=-0.015,cex.axis=cex.tickmark,las=1)  # y-axis
    axis(side=3,at=x.breaks+gap.pair[2],labels=x.breaks,tck=-0.015,cex.axis=cex.tickmark) # x-axis
    
    if(if.plot.legned){
      legend("topright",legend=labels[2],text.col=colors[2],cex=cex.legend,inset=0,box.col=FALSE) # legend
    }
    mtext(line=2,side=2,text=ylab,cex=cex.lab) # y-lab
    mtext(line=1.50,side=1,text=xlab,cex=cex.lab) # x-lab
    
  } else if (plot.method=="overlay"){
    ## overlay subplots
    if (if.plot.new) {
      par(mfrow=c(1,1),mar=c(1.5,4,1.8,0),omi=c(0.5,0.25,0.75,0.5),mgp=c(0,0.2,0))
    }
    ##1
    plot(x=1:length(x$v1)+gap.pair[1],
         y=x$v1,
         type=plot.type,
         ylim=ylim,
         xlim=xlim,
         ylab="",
         yaxt="n",
         xlab="",
         xaxt="n",
         col=colors[1],
         lwd=lwd.line
         )
    x2 <- 1:length(x$v2)+gap.pair[2]
    x2 <- x2 + (xlim[2]-xlim[1])*x.jitter
    points(x=x2,
           y=x$v2,
           type=plot.type,  
           col=colors[2],
           lwd=lwd.line
           )
    
    if(if.plot.legned)
      legend("topright",legend=labels,text.col=colors,cex=cex.legend,inset=0,box.col=FALSE)
    if (if.plot.gap){
      mtext(line=-.25,side=1,at=seq(xlim[1],xlim[2]),s.align,cex=3)
    }
    mtext(line=2.75,side=3,text=main,outer=FALSE,cex=2,adj=0)
    mtext(line=1.25,side=3,text=s.sub,outer=FALSE,cex=1.5,adj=0)
    mtext(line=2,side=2,text=ylab,cex=cex.lab)
    mtext(line=1.50,side=1,text=xlab,cex=cex.lab)

    axis(side=2,at=y.breaks,labels=y.breaks,tck=-0.006,cex.axis=cex.tickmark,las=1)
    axis(side=3,at=x.breaks+gap.pair[1],labels=x.breaks,col.axis=colors[1],tck=-0.006,cex.axis=cex.tickmark)
    axis(side=1,at=x.breaks+gap.pair[2],labels=x.breaks,col.axis=colors[2],tck=-0.006,cex.axis=cex.tickmark)
    
  }
  
  invisible()
  
}




## save.gmd <-
##   function(x,
##            pdfFpath="GMD_output.pdf",
##            pdfWidth=10,
##            pdfHeight=8,
##            if.save.text=TRUE,
##            ...)
## {
##   pdf(file=pdfFpath,width=pdfWidth,height=pdfHeight)
##   plot.gmd(x, if.plot.new=FALSE, ...)
##   dev.off()
##   print(sprintf("File is saved at: %s",file_path_as_absolute(pdfFpath)))
## }

