#' Recalculate locus metrics when individuals or populations are deleted from a genlight \{adegenet\} object
#'
#' When individuals are deleted from a genlight object generated by DArT, the locus metrics no longer
#' apply. For example, the Call Rate may be different considering the subset of individuals, compared
#' with the full set. This script recalculates those affected locus metrics, namely, avgPIC, CallRate,
#' freqHets, freqHomRef, freqHomSnp, OneRatioRef, OneRatioSnp, PICRef and PICSnp. Metrics that remain
#' unaltered are RepAvg and TrimmedSeq as they are unaffected by the removal of individuals.
#' 
#' The script optionally removes resultant monomorphic loci or loci
#' with all values missing and deletes them (using gl.filter.monomorphs.r). 
#' 
#' The script returns a genlight object with the recalculated locus metadata.
#'
#' @param x -- name of the genlight object containing SNP genotypes [required]
#' @param v -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return A genlight object with the recalculated locus metadata
#' @export
#' @author Arthur Georges (glbugs@@aerg.canberra.edu.au)
#' @examples
#' \dontrun{
#'    gl <- gl.recalc.metrics(testset.gl, v=2)
#' }
#' @seealso \code{\link{gl.filter.monomorphs}}


gl.recalc.metrics <- function(x, v=2){
  
  if (v > 0) {
    cat("Starting gl.recalc.metrics: Recalculating locus metrics\n")
  }

# Recalculate statistics
  x <- utils.recalc.avgpic(x,v=v)
  x <- utils.recalc.callrate(x,v=v)
  x <- utils.recalc.freqhets(x,v=v)
  x <- utils.recalc.freqhomref(x,v=v)
  x <- utils.recalc.freqhomsnp(x,v=v)
  
  if (v > 1) {  
    cat("Note: Locus metrics recalculated\n")
  }
  if (v > 0) {
    cat("Completed gl.recalc.metrics\n\n")
  }
  
  return <- x

}  