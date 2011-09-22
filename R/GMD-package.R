##' Computate Generalized Minimum Distance (GMD) between
##' discrete distributions
##'
##' \tabular{ll}{
##' Package: \tab GMD \cr
##' Type: \tab Package \cr
##' Version: \tab 0.2 \cr
##' Date: \tab Thu Sep 22 2011 \cr
##' License: \tab GPL (>= 2) \cr
##' }
##' This package contains functions for GMD computation, with
##' GMD algorithm implemented in C to interface with R.
##' 
##' To install from online repositories (e.g. CRAN),
##'
##' \code{install.packages(pkgs="GMD", repos="http://cran.r-project.org")}
##' 
##' To install from a downloaded source file,
##'
##' \code{install.packages(pkgs="GMD_<current-version>.tar.gz", repos=NULL)}
##' 
##' For a complete list of functions, use
##'
##' \code{library(GMD); ls("package:GMD")}
##'
##' @name GMD-package
##' @aliases GMD GMD-package
##' @docType package
##' @title The Package for Generalized Minimum Distance (GMD) Computation
##' @references Zhao et al (2011),
##' "Systematic Clustering of Transcription Start Site Landscapes",
##' \emph{PLoS ONE} \bold{6}(8): e23409.
##' \url{http://dx.plos.org/10.1371/journal.pone.0023409}
##'
##' See \code{citation("GMD")} for BibTeX entries for LaTeX users.
##' @author Xiaobei Zhao and Albin Sandelin
##' 
##' Maintainer: Xiaobei Zhao \email{xiaobei (at) binf.ku.dk}
##' @concept GMD gmd distance nonparametric optimize cluster classif
##' @keywords package
##' @seealso \code{\link{gmd}}, \code{\link{gmdm}}, \code{\link{cage}}
##' @examples
##' require(GMD) # load GMD
##' help(GMD) # a help document of GMD 
##' data(package="GMD") # a list of datasets available in GMD
##' ls("package:GMD") # a list of functions available in GMD
##' citation("GMD") # for citation
##' demo("GMD-demo") # run the demo
NULL

