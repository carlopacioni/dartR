#' Report summary of Read Depth for each locus
#'
#' SNP datasets generated by DArT report AvgCountRef and AvgCountSnp as counts of sequence tags for the reference and alternate alleles respectively.
#' These can be used to backcalculate Read Depth. Fragment presence/absence datasets as provided by DArT (SilicoDArT) provide Average Read Depth and 
#' Standard Deviation of Read Depth as stanard columns in their report.
#' 
#' Filtering on Read Depth using the companion script gl.filter.rdepth can be on the basis of loci with exceptionally low counts, 
#' or loci with exceptionally high counts.
#' 
#' The minimum, maximum and mean read depth are provided. Output also is a histogram of read depth, accompanied by a box and 
#' whisker plot presented either in standard (boxplot="standard") or adjusted for skewness (boxplot=adjusted). 
#' 
#' Refer to Tukey (1977, Exploratory Data Analysis. Addison-Wesley) for standard
#' Box and Whisker Plots and Hubert & Vandervieren (2008), An Adjusted Boxplot for Skewed
#' Distributions, Computational Statistics & Data Analysis 52:5186-5201) for adjusted
#' Box and Whisker Plots.
#' 
#' @param x -- name of the genlight object containing the SNP data [required]
#' @param boxplot -- if 'standard', plots a standard box and whisker plot; if 'adjusted',
#' plots a boxplot adjusted for skewed distributions [default 'adjusted']
#' @param range -- specifies the range for delimiting outliers [default = 1.5 interquartile ranges]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return -- dataframe with loci that are outliers
#' @importFrom graphics hist
#' @importFrom robustbase adjbox
#' @export
#' @author Arthur Georges (Post to \url{https://groups.google.com/d/forum/dartr})

# Last amended 3-Feb-19

gl.report.rdepth <- function(x, boxplot="adjusted", range=1.5, verbose=2) {

  # TIDY UP FILE SPECS
  
  build ='Jacob'
  funname <- match.call()[[1]]
  # Note does not draw upon or modify the loc.metrics.flags
  
  # FLAG SCRIPT START
  
  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }
  
  cat("Starting",funname,"[ Build =",build,"]\n")

  # STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    stop("Fatal Error: genlight object required!\n")
  }
  
  if (all(x@ploidy == 1)){
    cat("  Processing Presence/Absence (SilicoDArT) data\n")
  } else if (all(x@ploidy == 2)){
    cat("  Processing a SNP dataset\n")
  } else {
    stop("Fatal Error: Ploidy must be universally 1 (fragment P/A data) or 2 (SNP data)!")
  }
  
# DO THE JOB

  if (all(x@ploidy == 1)){
    rdepth <- x@other$loc.metrics$AvgReadDepth
  } else if (all(x@ploidy == 2)){
    rdepth <- x@other$loc.metrics$rdepth
  } 
  lower <- round(min(rdepth)-1,0)
  upper <- round(max(rdepth)+1,0)
  
  cat("No. of loci =", nLoc(x), "\n")
  cat("No. of individuals =", nInd(x), "\n")
  cat("  Miniumum read depth: ",round(min(rdepth),2),"\n")
  cat("  Maximum read depth: ",round(max(rdepth),2),"\n")
  cat("  Mean read depth: ",round(mean(rdepth),3),"\n\n")

  # Determine the loss of loci for a given filter cut-off
  retained <- array(NA,21)
  pc.retained <- array(NA,21)
  filtered <- array(NA,21)
  pc.filtered <- array(NA,21)
  percentile <- array(NA,21)
  for (index in 1:21) {
    i <- (index - 1)/20
    i <- (i - 1)*(1-upper) + 1
    percentile[index] <- i
    retained[index] <- length(rdepth[rdepth >= percentile[index]])
    pc.retained[index] <- round(retained[index]*100/nLoc(x),1)
    filtered[index] <- nLoc(x) - retained[index]
    pc.filtered[index] <- 100 - pc.retained[index]
  }
  df <- cbind(percentile,retained,pc.retained,filtered,pc.filtered)
  df <- data.frame(df)
  colnames(df) <- c("Threshold", "Retained", "Percent", "Filtered", "Percent")
  df <- df[order(-df$Threshold),]
  rownames(df) <- NULL
  #print(df)
  
  # Plot a histogram of read depth
  par(mfrow = c(2, 1),pty="m")
  
  # Prepare for plotting
  if (all(x@ploidy==2)){
    title <- paste0("SNP data (DArTSeq)\nRead Depth by locus")
  } else {
    title <- paste0("Fragment P/A data (SilicoDArT)\nRead Depth by locus")
  }  
  # Save the prior settings for mfrow, oma, mai and pty, and reassign
  op <- par(mfrow = c(2, 1), oma=c(1,1,1,1), mai=c(0.5,0.5,0.5,0.5),pty="m")
  # Set margins for first plot
  par(mai=c(1,0.5,0.5,0.5))
  # Plot Box-Whisker plot
  if (boxplot == "standard"){
    whisker <- boxplot(rdepth, horizontal=TRUE, col='red', range=range, main = title)
    if (length(whisker$out)==0){
      cat("  Standard boxplot, no adjustment for skewness\n")
      outliers <- NULL
    } else {
      outliers <- data.frame(Locus=as.character(x$loc.names[rdepth %in% whisker$out]),
                             ReadDepth=whisker$out
      )
      cat("  Standard boxplot, no adjustment for skewness\n")
    }
    
  } else {
    whisker <- robustbase::adjbox(rdepth,
                                  horizontal = TRUE,
                                  col='red',
                                  range=range,
                                  main = title)
    if (length(whisker$out)==0){
      cat("  Boxplot adjusted to account for skewness\n")
      outliers <- NULL
    } else {
      outliers <- data.frame(Locus=as.character(x$loc.names[rdepth %in% whisker$out]),
                             ReadDepth=whisker$out
      )
      cat("  Boxplot adjusted to account for skewness\n")
    }
  }  
  # Set margins for second plot
  par(mai=c(0.5,0.5,0,0.5))  
  # Plot Histogram
  hist(rdepth, col='red',breaks=100, main=NULL)
  
 # Output the outlier loci 
  if (length(whisker$out)==0){
    cat("  No outliers detected\n")
  } else {  
    if (verbose >=3){
      cat("  Outliers detected -- \n")
      print(outliers)
    }  
  }  
  
  # FLAG SCRIPT END
  
  if(verbose>=1){cat("Completed:",funname,"\n")}
  
  # Reset the par options    
  par(op)

  return(outliers)

}
