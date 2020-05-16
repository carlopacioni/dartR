#' Concatenates DArT trimmed sequences and outputs a nexus file (for BEAST).
#'
#' Concatenated sequence tags are useful for phylogenetic methods where
#' information on base frequencies and transition and transversion ratios are
#' required (for example, Maximum Liklihood methods). Where relevant,
#' heterozygous loci are resolved before concatenation by random allele
#' assignment.
#'
#' This function generates and aligment for use in BEAST or similar software
#' that accept a nexus file with sequence data. In these analyses, it is often
#' preferible to limit the analysis to reads that have several SNPs. For this
#' reason, the argument \code{min_nSNPs} can be used to retain only loci with a
#' minimum number of SNPs. Loci with multiple SNPs (i.e. secondaries) are
#' resolved in one sequence only.
#'
#' Currently only method 1 is employed
#'
#' Method 1 -- the heterozyous state is resolved by randomly assigning one or
#' the other SNP variant to the individual. The resultant sequence fragments are
#' concatenated across loci to generate a single composite haplotype to be used
#' in subsequent phylogenetic analyses.
#'
#' Trimmed sequences for which the SNP has been trimmed out, rarely, by adaptor
#' mis-identity are deleted.
#'
#' A charset block is appended at the end of the nexus file to facilitate
#' partitioning of the sequences.
#'
#' The script writes out the composite haplotypes for each individual as a nexus
#' file. Requires 'TrimmedSequence' to be among the locus metrics
#' (\code{@other$loc.metrics}) and information of the type of alleles (slot
#' loc.all e.g. "G/A") and the position of the SNP in slot position of the
#' ```genlight``` object (see testset.gl@position and testset.gl@loc.all for how
#' to format these slots.)
#'
#' @param gl -- name of the DArT genlight object [required]
#' @param method -- 1 | 2. Type method=0 for a list of options. Only mathod=1 implemented at the moment [method=1]
#' @param outfile -- name of the output file (fasta format) [output.fasta]
#' @param outpath -- path where to save the output file (set to tempdir by
#'   default)
#' @param min_nSNPs Minimum number of SNPs to retain a locus [min_nSNPs=3]
#' @param probar -- if TRUE, a progress bar will be displayed for long loops
#'   [default = TRUE]
#' @return A new gl object with all loci rendered homozygous
#' @export
#' @importFrom ape write.nexus.data
#' @importFrom utils combn edit flush.console getTxtProgressBar read.csv
#'   setTxtProgressBar txtProgressBar write.csv write.table
#' @importFrom methods new
#' @importFrom  stringr str_sub
#' @import data.table
#' @author Carlo Pacioni (Post to \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' gl <- gl.filter.repavg(testset.gl,t=1)
#' gl <- gl.filter.callrate(testset.gl,t=.98)
#' glBEAST <- gl2nexusBEAST(gl)


gl2nexusCompact <- function(gl, method=1, outfile="output.nex", outpath=tempdir(), 
                          min_nSNPs=3, probar=TRUE) {
  outfile <- file.path(outpath, outfile)
  
  if(!is(gl, "genlight")) {
    stop("Fatal Error: Specify a genlight object\n")
  }
  if(length(gl@other$loc.metrics$TrimmedSequence) != nLoc(gl)) {
    stop("Fatal Error: Data must include Trimmed Sequences for each loci in a column called 'TrimmedSequence' in the @other$loc.metrics slot.\n")
  }
  if(length(gl@position) != nLoc(gl)) {
     stop("Fatal Error: Data must include position information for each loci in the @position slot.\n")
  }
  if(length(gl@loc.all) != nLoc(gl)) {
    stop("Fatal Error: Data must include type of alleles in the @loc.all slot.\n")
  }
  
  if(method==1){
    cat("Randomly selecting one allele in heterozygotes (1), concatenating trimmed sequence\n")
  #} else if (method==3) {
  #  cat("Assigning ambiguity codes to heterozygote SNPs, concatenating SNPs\n")
  #}else if (method==4) {
  #  cat("Randomly allocating heterozygotes (1) to homozygotic state (0 or 2), concatenating SNPs\n")
  } else {
    cat("Fatal Error: Parameter method out of range.\n")
    cat("Replace score for heterozygotic loci with"):
    cat("  method=1 -- Randomly selecting one allele in heterozygotes (1), concatenating trimmed sequence [default]\n")
    stop()
  }

# METHOD = RANDOM ASSIGNMENT

if (method==1) {

# Randomly allocate heterozygotes (1)
  matrix <- t(as.matrix(gl))
  #cat("Randomly allocating heterozygotes (1) to homozygote state (0 or 2)\n")
  #pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  #getTxtProgressBar(pb)
  
  sample_nms <- colnames(matrix)
  locMetr <- data.table(gl@other$loc.metrics, keep.rownames = TRUE)
  locMetr[, lenTrim := nchar(as.character(TrimmedSequence))]
  index <- locMetr[,lenTrim] > gl@position
  if (sum(index)!=nLoc(gl)) {
    message(paste("Not all SNP positions are within the length of the trimmed sequences. Those loci will be deleted (",sum(!index),").\n"  ) )
    message("Details of loci being deleted\n")
    message(paste(locMetr[!index, uid], "\n"))
  }
  locMetr <- locMetr[index, ]
  matrix <- matrix[index, ]
  
  nSNPs <- locMetr[, .N, by=clone]
 
  locMetr <- merge.data.table(locMetr, nSNPs, by="clone")
  matrix <- matrix[locMetr[, N>=min_nSNPs],]
  locMetr <- locMetr[N>=min_nSNPs,]
  if(nrow(matrix) == 0) {
    stop(paste("No loci met the threshold min_nSNPs:", min_nSNPs))
  }
  
  r <- nrow(matrix)
  c <- ncol(matrix)
  system.time(
  for (i in 1:r) {
    for (j in 1:c) {
      if (matrix[i,j] == 1 && !is.na(matrix[i,j])) {
        # Score it 0 or 2
        matrix[i,j] <- sample(c(0, 2), 1)
      }
    }
  #  setTxtProgressBar(pb, i/r)
  }
  )
 
  
# Prepare the output fastA file
  cat("Generating haplotypes ... This may take some time\n")

#  sink(outfile)

# For each individual, and for each locus, generate the relevant haplotype 
#  seq <- rep(" ", c)
  #pb <- txtProgressBar(min=0, max=1, style=3, initial=0, label="Working ....")
  #getTxtProgressBar(pb)
  
  # Shift the index for snppos to start from 1 not zero only once
  locMetr[, SnpPosition := SnpPosition + 1]
  
  clones <- locMetr[, unique(clone)]
  message(paste("Number of loci >=", min_nSNPs, ":", length(clones)))
  
  ncln <- 0
  mseqs <- matrix(" ", nrow=length(clones), ncol=c) 
  system.time(
  for(cln in clones) {
    ncln <- ncln + 1
    m_temp <- matrix[locMetr[, clone] == cln,]
    locMetr_temp <- locMetr[clone == cln, ]
    
    for (j in seq_len(ncol(m_temp))) {
      trimmed <- as.character(locMetr_temp[1, TrimmedSequence])
      # If the SNP state is missing in all positions, assign NNNNs  
      if (sum(is.na(m_temp[,j])) == nrow(m_temp)) {
        mseqs[ncln, j] <- str_pad("N", nchar(trimmed), side = c("right"), pad="N")
          next
      }
      for (i in seq_len(nrow(m_temp))) {
        snpos <- locMetr_temp[i, SnpPosition]
          # If the score is homozygous for the reference allele
          if (m_temp[i,j] == 0 && !is.na(m_temp[i,j])) { 
            next
          } else {
            if (m_temp[i,j] == 2 && !is.na(m_temp[i,j])) {
              snpbase <- str_sub(rownames(m_temp)[i], start = -1, end = -1) # this could be done with substr without the need for an additional package
              
              # Change the SNP state to the alternate
              str_sub(trimmed, start=(snpos), end=(snpos)) <- snpbase
            } else {
              if (is.na(m_temp[i,j])) {
                str_sub(trimmed, start=(snpos), end=(snpos)) <- "N"
              } else {
                stop(paste("Something went wrong allele values for locus", 
                           rownames(m_temp)[i], "and sample", colnames(m_temp)[j],
                           "is not 0, 1 or 2"))
              }
            }
          }
      }  
      mseqs[ncln, j] <- trimmed
    }
  }
  )
      #setTxtProgressBar(pb, i/r)
      seqs_vec <- apply(X = mseqs, MARGIN = 2, FUN = paste, collapse="")
      names(seqs_vec) <- sample_nms
      ape::write.nexus.data(x = seqs_vec, file = outfile)
      
      # write charset block
      pos <- 0
      charsets <- character(length(clones))
      
      for(cln in clones) {
        charsets[which(clones == cln)] <- paste("charset", cln, "=", 1 + pos, "-", 
                             pos + locMetr[clone == cln, lenTrim][1], ";")
        pos <-  pos + locMetr[clone == cln, lenTrim][1]
      }
      
      write(c("\n","begin assumptions;", charsets, "end;"), 
            file = outfile, append = TRUE, sep = "\n")
      
}
glnew <- new("genlight", gen=t(matrix), ind.names=sample_nms, 
             loc.names=rownames(matrix), 
             loc.all=str_sub(rownames(matrix), start=nchar(rownames(matrix)) - 2, end=-1),
             position=locMetr[, SnpPosition - 1], other=list(loc.metrics=locMetr))
return(glnew)
}


