\name{chipseq}
\title{ChIP-seq data: ChIP-seq Enrichment around TSSs}
\usage{chipseq_mES
chipseq_hCD4T}
\description{The Distributions of Histone Modification Enrichment (and Others) by ChIP-seq reads that are binned,
aligned and averaged around +/-5000nt of Transcription Start Sites (TSSs) of scattered-type TSSDs
(see References).}
\details{\code{chipseq_mES} is a list of 6 named ChIP-seq read distributions from mouse ES cells.\cr
\code{chipseq_hCD4T} is a list of 40 named ChIP-seq read distributions from human CD4+ T cells.\cr}
\alias{chipseq}
\alias{chipseq_mES}
\alias{chipseq_hCD4T}
\docType{data}
\references{Zhao et al (2011),
"Systematic Clustering of Transcription Start Site Landscapes",
\emph{PLoS ONE} \bold{6}(8): e23409.

\url{http://dx.plos.org/10.1371/journal.pone.0023409}}
\seealso{\code{\link{gmdp}} and \code{\link{gmdm}}, with examples using \code{\link{chipseq}}. \code{\link{cage}} for data of Transcription Start Sites (TSSs) by CAGE tags.}
\keyword{datasets}
\examples{require(GMD)
help(chipseq)
data(chipseq_mES)
class(chipseq_mES)
length(chipseq_mES)
names(chipseq_mES)
\dontrun{chipseq_mES}

data(chipseq_hCD4T)
names(chipseq_hCD4T)}
