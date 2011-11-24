\name{ts2df}
\alias{ts2df}
\title{Convert time series to data frame}
\usage{ts2df(x, mon1=1, addYr=FALSE, omit=FALSE)}
\description{A copy of wq:::ts2df; see \code{ts2df} in package:wq for details}
\details{see \code{ts2df} in package:wq for details.
Note: wq_0.3-4 asks for R (>= 2.12.0); but GMD supports R (>= 2.9.0).}
\arguments{\item{x}{monthly time series vector}
\item{mon1}{starting month number, i.e., first column of the data frame}
\item{addYr}{rows are normally labelled with the year of the starting
month, but \code{addYr = TRUE} will add 1 to this year number}
\item{omit}{if \code{TRUE}, then rows with any \code{NA} will be removed.}}
