library(rtracklayer)

GRangeCreation <- function(bamfile, chname = T){
  gr_bamfile <- readGAlignments(bamfile, use.names = T)
  if(chname){
    seqlevels(gr_bamfile) <- sub('chr', '', seqlevels(gr_bamfile)) #remove 'chr' from contig names
  }
  
  return(gr_bamfile)
} 

summarize <- function(signal, n, FUN=mean) {
  N <- length(signal)
  if( N == n ) return(signal)
  if( N < n ) {
    eachvals <- rep(floor(n/N), N)
    residual <- n - sum(eachvals)
    if( residual > 0 ) eachvals[1:residual] <- eachvals[1:residual] + 1
    return(unlist(lapply(1:length(eachvals), function(i) rep(signal[i], eachvals[i]))))
  } else {
    eachvals <- rep(floor(N/n), n)
    residual <- N - sum(eachvals)
    if( residual > 0 ) eachvals[1:residual] <- eachvals[1:residual] + 1
    splitvals <- unlist(lapply(1:length(eachvals), function(i) rep(i, eachvals[i])))
    return(tapply(signal, splitvals, FUN, na.rm=TRUE))
  }
}

coverageFromFile <- function(txdb, file, control_file = NULL, gr, Nbin=NULL, ignore.strand=FALSE, invertBamStrand=FALSE, revMinusStrand=TRUE, strand, change_chr_names = T) {
  
  if( grepl('\\.bam$', tolower(file[1])) ) {
    gr_bamfile <- GRangeCreation(file, chname = change_chr_names)
    strand(gr_bamfile) = strand
    #gr_bamfile = subsetByOverlaps(gr_bamfile, exons(txdb), invert = T)
    sampl_cov = coverage(gr_bamfile)
    coverage_unstr = sampl_cov
  }
  else if( grepl('\\.bw$', tolower(file[1])) ) {
    
  }
  covRLE <- coverageFromRLE(coverage_unstr, gr, Nbin)
  return(covRLE)

}

checkCompatibleSeq <- function(coverage, gr) {
  bigWigSeqLengths <- sapply(coverage, length)
  annotationSeqLengths <- seqlengths(seqinfo(gr))
  commonSeq <- intersect(names(bigWigSeqLengths), names(annotationSeqLengths))
  if( length(commonSeq) == 0 ) 
    stop('coverage and gr have no common sequences')
  #if( !identical(bigWigSeqLengths[commonSeq], annotationSeqLengths[commonSeq]) ) 
  #  stop('coverage and gr file have different sequence lengths')
  return(commonSeq)
}

coverageFromRLE <- function(Rle_list, gr, Nbin=NULL) {
  
  grCoverage <- lapply(checkCompatibleSeq(Rle_list, gr), function(chr) {
    gr_chr <- gr[seqnames(gr) == chr]
    chr_cov <- Views(Rle_list[[chr]], ranges(gr_chr))
    lapply(chr_cov, as.numeric)
  })
  names(grCoverage) <- NULL
  grCoverage <- do.call('c', grCoverage)
  # give the initial ordering
  grCoverage <- grCoverage[names(gr)]
  if( !is.null(Nbin) ) grCoverage <- lapply(grCoverage, summarize, n=Nbin)
  return(grCoverage)
  
}

coverageFromBW <- function(bwfile, gr, Nbin=NULL) {
  coverage <- import(bwfile, as = 'RleList',which=gr)
  # check compatibility between BigWig and annotation
  bigWigSeqLengths <- sapply(coverage, length)
  annotationSeqLengths <- seqlengths(seqinfo(gr))
  commonSeq <- intersect(names(bigWigSeqLengths), names(annotationSeqLengths))
  if( length(commonSeq) == 0 ) 
    stop('quantifyExpressionsFromBWs: annotation and bigwig file have no common sequences')
  
  grCoverage <- lapply(names(coverage), function(chr) {
    gr_chr <- gr[seqnames(gr) == chr]
    chr_cov <- Views(coverage[[chr]], ranges(gr_chr))
    lapply(chr_cov, as.numeric)
  })
  names(grCoverage) <- NULL
  grCoverage <- do.call('c', grCoverage)
  # give the initial ordering
  grCoverage <- grCoverage[names(gr)]
  if( !is.null(Nbin) ) grCoverage <- lapply(grCoverage, summarize, n=Nbin)
  return(grCoverage)
}





