
## ------------------------------------------------------------------------
## Wrapping the Library Routine
## ------------------------------------------------------------------------
.mdpa <- function(v1, v2){
  if (length(v1) != length(v2)){
    stop("The `mdpa` function should take two equal-length vectors as input!\n")
  }
  if (sum(v1) != sum(v2)){
    stop("The `mdpa` function should take two equal-mass vectors as input!\n")
  }
  .Call("mdpa", v1, v2, PACKAGE="GMD")
}


.gmd0 <- function(v1, v2, pseudocount=0){
  if (length(v1) != length(v2)){
    stop("The `gmd0` function should take two equal-length vectors as input!\n")
  }
  res <- .Call("gmd0", v1, v2, pseudocount, PACKAGE="GMD")
  return(res)
}


.gmd <- function(v1, v2, pseudocount=0){
  res <- .Call("gmd", v1, v2, pseudocount, PACKAGE="GMD")
  return(res)
}


## ------------------------------------------------------------------------
## Other internal functions
## ------------------------------------------------------------------------

.wordwrap <-
  function(x,len)
{
  l <- nchar(x)
  m <- matrix(ncol=2,nrow=ceiling(l/len))
  for (i in 1:nrow(m)){
    m[i,] <- c(1+len*(i-1),min(len*i,l))
  }
  res <- mapply(FUN=substr,m[,1],m[,2],MoreArgs=list(x=x))
  paste(res,sep="",collapse="\n")
}

.is.gmd <-
  function(x)
{
  "gmd" %in% class(x)
}

.is.gmdm <-
  function(x)
{
  "gmdm" %in% class(x)
}


## ------------------------------------------------------------------------
## Other params
## ------------------------------------------------------------------------
.brewer.pal.Dark2 <-
  c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")


