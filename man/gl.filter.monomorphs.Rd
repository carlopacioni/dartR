% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl.filter.monomorphs.r
\name{gl.filter.monomorphs}
\alias{gl.filter.monomorphs}
\title{Remove monomorphic loci, including those with all NAs}
\usage{
gl.filter.monomorphs(x, verbose = 2, pb = FALSE)
}
\arguments{
\item{x}{-- name of the input genlight object [required]}

\item{verbose}{-- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]}

\item{pb}{-- display progress bar [FALSE]}
}
\value{
A genlight object with monomorphic loci removed
}
\description{
This script deletes monomorphic loci from a genlight \{adegenet\} object
}
\details{
A DArT dataset will not have monomorphic loci, but they can arise when populations are deleted by assignment or by using
the delete option in gl.pop.recode(). Retaining monomorphic loci unnecessarily increases the size of the dataset.
}
\examples{
gl <- gl.filter.monomorphs(testset.gl, verbose=3)
}
\author{
Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
}
