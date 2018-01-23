#' A utility script to recalculate the callrate by locus after some populations have been deleted
#'
#' SNP datasets generated by DArT have missing values primarily arising from failure to call a SNP because of a mutation
#' at one or both of the the restriction enzyme recognition sites. The locus metadata supplied by DArT has callrate included,
#' but the call rate will change when some individuals are removed from the dataset. This script recalculates the callrate
#' and places these recalculated values in the appropriate place in the genlight object.
#'
#' @param gl -- name of the genlight object containing the SNP data [required]
#' @param v -- v=0, silent; v=1, low verbosity; v=2, high verbosity [default 1]
#' @return The modified genlight object
#' @author Arthur Georges (glbugs@aerg.canberra.edu.au)
#' @examples
#' result <- utils.recalc.callrate(testset.gl)

utils.recalc.callrate <- function(gl, v=1) {
 x <- gl
   
  if(class(x) == "genlight") {
     #cat("Modifying a genlight object\n")
   } else {
     cat("Fatal Error: Specify a genlight object\n")
     stop()
  }

  # Do the deed
     x@other$loc.metrics$CallRate <- 1-(glNA(gl,alleleAsUnit=FALSE))/nInd(gl)

   if (v>0) {cat("Callrate recalculated\n")}
   
   return(x)
}

