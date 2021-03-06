% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.subsample.loci.r
\name{utils.subsample.loci}
\alias{utils.subsample.loci}
\title{Subsample n loci from a genlight object and return as a genlight object}
\usage{
utils.subsample.loci(x, k, method = "random", verbose = 2)
}
\arguments{
\item{x}{-- name of the genlight object containing the SNP genotypes by specimen and population [required]}

\item{k}{-- number of loci to include in the subsample [required]}

\item{method}{-- "random", in which case the loci are sampled at random; or avgPIC, in which case the top n loci
ranked on information content (AvgPIC) are chosen [default "random"]}

\item{verbose}{-- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]}
}
\value{
A genlight object with k loci
}
\description{
This is a support script, to subsample a genlight \{adegenet\} object based on loci. Two methods are used
to subsample, random and based on information content (avgPIC).
}
\examples{
#result <- dartR:::utils.subsample.loci(testset.gl, k=200, method="avgPIC")
}
\author{
Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})
}
