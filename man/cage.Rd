\name{cage}
\title{CAGE Data: Transcription Start Site Distributions (TSSDs) by CAGE tags}
\usage{cage
cagel}
\description{Transcription Start Site Distributions (TSSDs) by CAGE tags.}
\details{\code{cage} is a list of 20 named TSSDs.
\code{cagel} is a longer version of \code{\link{cage}}, with 50 named TSSDs.}
\alias{cage}
\alias{cagel}
\docType{data}
\references{Zhao et al (2011),
"Systematic Clustering of Transcription Start Site Landscapes",
\emph{PLoS ONE} \bold{6}(8): e23409.

\url{http://dx.plos.org/10.1371/journal.pone.0023409}}
\seealso{\code{\link{gmdp}} and \code{\link{gmdm}}, with examples using \code{\link{cage}}. \code{\link{chipseq}} for histone marks by ChIP-seq reads.}
\keyword{datasets}
\examples{help(cage)
data(cage)
class(cage)
length(cage)
names(cage)
\dontrun{cage}

data(cagel)
names(cagel)
\dontrun{cagel}}
