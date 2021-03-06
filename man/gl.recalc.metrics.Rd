% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.recalc.metrics.r
\name{gl.recalc.metrics}
\alias{gl.recalc.metrics}
\title{Recalculate locus metrics when individuals or populations are deleted from a genlight \{adegenet\} object}
\usage{
gl.recalc.metrics(x, verbose = 2)
}
\arguments{
\item{x}{-- name of the genlight object containing SNP genotypes [required]}

\item{verbose}{-- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]}
}
\value{
A genlight object with the recalculated locus metadata
}
\description{
When individuals are deleted from a genlight object generated by DArT, the locus metrics no longer
apply. For example, the Call Rate may be different considering the subset of individuals, compared
with the full set. This script recalculates those affected locus metrics, namely, avgPIC, CallRate,
freqHets, freqHomRef, freqHomSnp, OneRatioRef, OneRatioSnp, PICRef and PICSnp. Metrics that remain
unaltered are RepAvg and TrimmedSeq as they are unaffected by the removal of individuals.
}
\details{
The script optionally removes resultant monomorphic loci or loci
with all values missing and deletes them (using gl.filter.monomorphs.r).

The script returns a genlight object with the recalculated locus metadata.
}
\examples{
  gl <- gl.recalc.metrics(testset.gl, verbose=2)
}
\seealso{
\code{\link{gl.filter.monomorphs}}
}
\author{
Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
}
