library(rtracklayer)
h3k4 = import("/hpcnfs/data/GN2/gmandana/annotation/HeLa_H3K4me3_peaks_bed6.bed")
seqlevels(h3k4) = sub("chr", "", seqlevels(h3k4))
h3k4 = reduce(h3k4, min.gapwidth = 1000)

GRangeCreation <- function(bamfile, chname = T){
  gr_bamfile <- readGAlignments(bamfile, use.names = T)
  if(chname){
    seqlevels(gr_bamfile) <- sub('chr', '', seqlevels(gr_bamfile)) #remove 'chr' from contig names
  }
  
  return(gr_bamfile)
} 

txMatrixCreation <- function(annotation_switch, txdb, tx_filtered, tx_names, gr_bamfile, non_zero_counts, tss_upstream, tss_downstream, 
                             magnitude_factor, analysis, tss_proximal, tss_distal, tss, tss_analysis){
  ###warning
  
  if(sum(seqlevels(selectedtx) %in% seqlevels(txdb)) == 0){
    stop('selectedtx and txdb have no common sequences')
  }
  if(sum(seqlevels(selectedtx) %in% seqlevels(gr_bamfile)) == 0){
    stop('selectedtx and bamfile have no common sequences')
  }
  
  # Creation
  transcripts_matrix <- as.data.frame(matrix(FALSE, length(tx_filtered), 7,
                                             dimnames = list(tx_names, 
                                                             c('biotype','width','active', 'h3k4me3',
                                                               'no_other_active_tss','no_unann_tss','coverage'))))
  if(annotation_switch){
    transcripts_matrix$gene_name <- sapply(strsplit((names(tx_filtered)), '\\.'), function(x){x[1]})
    transcripts_matrix[,'biotype'] <- factor(AnnotationDbi::select(txdb, tx_names, 'TXTYPE', 'TXNAME')$TXTYPE)
  }else{
    transcripts_matrix$gene_name <- tx_filtered$index
    transcripts_matrix[,'biotype'] <- "noncoding"
  }
  # 2) Select transcripts based on their length (40kbp)
  transcripts_matrix[,'width'] <- width(tx_filtered) # > analysis_down_width
  
  # 3) I am putting a filtering step here : coverage based
  transcripts_matrix[,'coverage']  = countOverlaps(analysis, gr_bamfile)
  
  # 4) Identify active TSS (Select genes where the value 3′ from the TSS was at least 10 times higher than the value 5′ of the TSS)
  tss_upstream_counts <- countOverlaps(tss_upstream, gr_bamfile) # 1 is added to avoid problems in the ratio
  tss_downstream_counts <- countOverlaps(tss_downstream, gr_bamfile) # 1 is added to avoid problems in the ratio
  
  transcripts_matrix[,'active'] <- floor((non_zero_counts+tss_downstream_counts)/(non_zero_counts+tss_upstream_counts)) # > magnitude_factor
  
  transcripts_matrix[,'h3k4me3'] = countOverlaps(tss_analysis, h3k4)
  
  # 5) identify additional active TSS within analysis range
  
  unique_tss_active <- reduce(unique(tss[transcripts_matrix[,'active']>magnitude_factor]), min.gapwidth = 1000) # tsses within 1kb are merged
  transcripts_matrix[,'no_other_active_tss'] <- countOverlaps(analysis, unique_tss_active) # == 1
  
  # 6) rejected if the TSS-proximal signal was not more than magnitude factr times the distal signal
  tss_proximal_counts <- countOverlaps(tss_proximal, gr_bamfile)
  tss_distal_counts <- countOverlaps(tss_distal, gr_bamfile)
  transcripts_matrix[,'no_unann_tss'] <- floor((non_zero_counts+tss_proximal_counts)/(non_zero_counts+tss_distal_counts)) # >magnitude_factor
  
  return(transcripts_matrix)
}

FilteringOperation <- function(transcripts_matrix, ranges, minGeneLength = 40000, magnitude_factor = 10, coverage_factor = 150, sample, strand){
  
  
  transcripts_matrix_filtered1 = transcripts_matrix
  # 2) Select transcripts based on their length (40kbp/150kbp), here 40kbp
  transcripts_matrix_filtered2 <- subset(transcripts_matrix_filtered1, transcripts_matrix_filtered1$width >= minGeneLength)
  
  # 3)    COVERAGE  
  transcripts_matrix_filtered3 = subset(transcripts_matrix_filtered2, transcripts_matrix_filtered2$coverage >= coverage_factor)
  
  # 3)  Identify active TSS (Select genes where the value 3′ from the TSS was at least 10 times higher than the value 5′ of the TSS)
  transcripts_matrix_filtered4 <- subset(transcripts_matrix_filtered3, transcripts_matrix_filtered3$active >= magnitude_factor)
  
  # 4) identify additional active TSS within analysis range
  # I made the change to merge TSSes within a kb of each other
  
  #transcripts_matrix_filtered5 <- subset(transcripts_matrix_filtered4, transcripts_matrix_filtered4$h3k4me3 == 1 & transcripts_matrix_filtered4$no_other_active_tss == 1)
  transcripts_matrix_filtered5 <- subset(transcripts_matrix_filtered4, transcripts_matrix_filtered4$h3k4me3 == 1)
  
  # 5) rejected if the TSS-proximal signal was not more than 10 times the distal signal
  transcripts_matrix_filtered6 <- subset(transcripts_matrix_filtered5, transcripts_matrix_filtered5$no_unann_tss >= magnitude_factor)
  #transcripts_matrix_filtered6 = transcripts_matrix_filtered5
  
  transcripts_matrix_filtered6$txname = rownames(transcripts_matrix_filtered6)
  
  transcripts_matrix_filtered_final = as.data.frame(transcripts_matrix_filtered6 %>% group_by(gene_name) %>% dplyr::filter(row_number()==1)) # CHECK MY NOTEBOOK 
  rownames(transcripts_matrix_filtered_final) = transcripts_matrix_filtered_final$txname
  
  filteringstats <<- data.frame("initial number" = length(rownames(transcripts_matrix_filtered1)),
                                "length filter" = length(rownames(transcripts_matrix_filtered2)),
                                "coverage filter" = length(rownames(transcripts_matrix_filtered3)),
                                "active TSS filter" = length(rownames(transcripts_matrix_filtered4)),
                                "check for other TSS in analysis region" = length(rownames(transcripts_matrix_filtered5)),
                                "no unann TSS" = length(rownames(transcripts_matrix_filtered6)))
  
  #write.csv(filteringstats, file = paste(sample, strand, "_filtering.csv", sep = ""))
  
  return(transcripts_matrix_filtered_final)
}

# define ranges, binned bru coverage, remove and interpolate the exon signal
DataProcessing <- function(annotation_switch, bruTxNames, bamfile, control_bamfile, beforeTSS = 10000, geneLength = 150000, Nbins = 640, txdb, strand){
  # define ranges
  #bruTxNames <- BruDRBseq$Ensembl.Transcript.ID
  
  if(annotation_switch){
    tx_annotation <- transcripts(txdb, use.names=T)[bruTxNames]
    bruTxNames <- names(tx_annotation[width(tx_annotation)>geneLength])
    bruGR <- promoters(txdb, upstream = beforeTSS, downstream = geneLength)[bruTxNames]
    exonsBruGR <- reduce(exonsBy(txdb, 'tx', use.names=TRUE)[bruTxNames])
    #seqlevels(bruGR) <- seqlevels(exonsBruGR) <- paste0('chr',seqlevels(exonsBruGR))
  }else{
    bruGR <- promoters(bruTxNames, upstream = beforeTSS, downstream = geneLength)
  }
  # binned bru coverage 
  bruCovBin <- do.call('rbind', coverageFromFile(txdb, bamfile, NULL, bruGR, Nbin = Nbins, ignore.strand=T, invertBamStrand=FALSE, revMinusStrand=F, strand))
  
  #control file coverage
  bruCovBin_ctl <- do.call('rbind', coverageFromFile(txdb, control_bamfile, NULL, bruGR, Nbin = Nbins, ignore.strand=T, invertBamStrand=FALSE, revMinusStrand=F, strand))

    # remove and interpolate the exon signal if annotation present
  if(annotation_switch){
    
    bruGR.isexon.Nbin <- t(sapply(GRbaseIsExon(bruGR, exonsBruGR), summarize, n=Nbins))
    bruCovBin.NA <- bruCovBin
    bruCovBin.NA[bruGR.isexon.Nbin > 0] <- NA
    bruCovBin.imputed <- t(sapply(1:nrow(bruCovBin.NA), function(i) impute_na_tc(1:ncol(bruCovBin.NA), bruCovBin.NA[i,])))
    rownames(bruCovBin.imputed) <- rownames(bruCovBin)
    bruCovBin.imputed = bruCovBin
    
  }else{
    #assign, just make a copy, don't change the obj cos maybe the code breaks later lol
    bruCovBin.imputed = bruCovBin
  }
  #test_coverage <<- bruCovBin
  lst <- list(
    bruCovBin = bruCovBin,
    bruCovBin_ctl = bruCovBin_ctl,
    ranges = bruGR # the actual ranges being considered
  )
  return(lst)
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

impute_na_tc <- function(tpts, tcdata) {
  
  # impute NA values in a time course data as a linear model between non-NA values
  tc_impute_NA_linearmodel <- function(tpts, tcdata) {
    for( j in seq_along(tcdata) ) {
      if( is.na(tcdata[j]) ) {
        lower_boundary_j <- j-1
        higher_boundary_j <- j+1
        while( is.na(tcdata[higher_boundary_j]) & higher_boundary_j <= length(tcdata) ) {
          higher_boundary_j <- higher_boundary_j + 1
        }
        if( lower_boundary_j > 0 & higher_boundary_j <= length(tcdata) ) 
          if ( is.finite(tcdata[lower_boundary_j]) ) {
            x <- tpts[c(lower_boundary_j,higher_boundary_j)]
            y <- tcdata[c(lower_boundary_j,higher_boundary_j)]
            tcdata[(lower_boundary_j+1):(higher_boundary_j-1)] <- predict(lm(y ~ x), 
                                                                          data.frame(x=tpts[(lower_boundary_j+1):(higher_boundary_j-1)]))
          }
      }
    }
    return(tcdata)
  }
  
  # impute NA values in a time course from boundary values
  tc_impute_NA_boundaries <- function(tcdata) {
    
    forward_direction <- function(tcdata) {
      if( is.na(tcdata[1]) & !all(is.na(tcdata)) ) {
        lower_boundary_j <- higher_boundary_j <- 1
        while( is.na(tcdata[higher_boundary_j] & higher_boundary_j < length(tcdata) ) ) {
          higher_boundary_j <- higher_boundary_j + 1
        }
        tcdata[lower_boundary_j:(higher_boundary_j-1)] <- tcdata[higher_boundary_j]
      }
      return(tcdata)
    }
    
    tcdata <- forward_direction(tcdata)
    tcdata <- rev(forward_direction(rev(tcdata)))
    return(tcdata)
  }
  
  tcdata <- tc_impute_NA_linearmodel(tpts, tcdata)
  tcdata <- tc_impute_NA_boundaries(tcdata)
  return(tcdata)
  
}

GRbaseIsExon <- function(gr, exons_ann, ignore.strand=FALSE) {
  gr.isexon <- lapply(1:length(gr), function(i) {
    r1 <- IRanges(start=start(gr[i]):end(gr[i]), width = 1) # makes a ranges obj for each position in the transcript
    isexon <- countOverlaps(r1, ranges(exons_ann[[i]]))     # counts overlap with each exon
  })
  if( !ignore.strand ) {
    gr.isexon <- lapply(1:length(gr.isexon), function(i) {
      if( as.character(strand(gr[i])) == '+' ) return(gr.isexon[[i]]) else rev(gr.isexon[[i]])
    })
  }
  return(gr.isexon)
}

postpr = function(obj){
  
  #want to remove nonsense
  indi = obj$res$start_2 > -30 & !is.infinite(obj$res$start_2)
  obj$res = obj$res[indi,]; obj$ranges = obj$ranges[indi]; obj$coverage = obj$coverage = obj$coverage[indi,]; obj$fitted = obj$fitted[indi,]
  
  obj$res$quartile = cut(obj$res$elong, quantile(obj$res$elong, probs = seq(0,1,0.25)), labels = F, include.lowest = T)
  
  #get the sequence for gc
  adjranges = adj(obj$ranges, obj$res$start_2, kbperbin, analysis_up_width)
  first1k = getSeq(genome, promoters(adjranges, upstream = 0, downstream = 1000))
  obj$res$GCinit = rowSums(alphabetFrequency(first1k)[,c(2,3)]/width(first1k)*100)
  obj$res$GCquart = cut(obj$res$GCinit, quantile(obj$res$GCinit, probs = seq(0,1,0.25)), labels = F, include.lowest = T)
  
  return(obj)
  
}

nainorm = function(matrix){
  maxvals = unlist(lapply(1:nrow(matrix), function(x){max(matrix[x,])}))
  return(matrix/maxvals)
}

eststop = function(matrix, thr = 0.1, start = 50){
  #estimate when the signal drops below thresh; pass normaslised coverage matrix
  #assume signal spike start ~40 bins in, so start measurement at 40th bin -so add 40 to the final index
  
  #for ratio, I think the thresh is 1.0
  #this relies on a good filtering of the transcripts
  
  stops = lapply(1:nrow(matrix), function(x){which(as.vector(unlist(slide(matrix[x,][start:length(matrix[x,])], mean, .before = 5, .after = 5))) < thr)[1]})
  
  return(unlist(stops) + start)
}

reverse_cov = function(obj){
  
  test = do.call('rbind', lapply(1:nrow(obj$sample_cov), function(x){rev(obj$sample_cov[x,])}))
  rownames(test) = rownames(obj$sample_cov)
  test2 = do.call('rbind', lapply(1:nrow(obj$control_cov), function(x){rev(obj$control_cov[x,])}))
  rownames(test2) = rownames(obj$control_cov)
  obj$sample_cov = test
  obj$control_cov = test2
  return(obj)
  
}

evenmoreplots = function(filepath, obj, sample, analysis_up_width, analysis_down_width){
  
  pdf(filepath)
  
  kek = pheatmap(obj$coverage, cluster_cols = F, cutree_rows = 4, clustering_method="ward.D2", show_rownames = F,
                 clustering_distance_rows = dist(obj$coverage[,41:200]),
                 labels_col = -((analysis_up_width)/1000):((analysis_down_width)/1000), main = paste("sample", sample, sep =" "))
  
  hcl = hclust(dist(obj$coverage), method="complete")
  clu4 = cutree(hcl, k = 4)
  clusters = as.data.frame(clu4)
  clusters$ids = rownames(clusters)
  
  iwant = c("protein_coding", "lincRNA", "antisense")
  hah = obj$res[obj$res$biotype %in% iwant,]
  
  
  caption = paste("counts = ", lapply(hah %>% dplyr::count(biotype) %>% dplyr::select(n), toString), sep = "")
  
  one = ggplot(hah, aes(x = factor(biotype), y = elong)) + 
    geom_boxplot(aes(color=factor(biotype)), width = 0.5, lwd=1.0, notch=T) + 
    labs(x = "Biotype", y = "Rate of Elongation", title = "Elongation Rate by Biotype", caption = caption) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylim(0,3000)
  
  print(one)
  
  three = ggplot(hah, aes(x = factor(biotype), y = elong)) + 
    geom_boxplot(aes(color=factor(biotype)), width = 0.5, lwd=1.0, notch=T) + 
    facet_grid(~GCquart) +
    labs(x = "Biotype", y = "Rate of Elongation", title = "Elongation Rate, quartiles of GC") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ylim(0,3000)
  
  print(three)
  
  four = ggplot(hah, aes(x = GCinit, y = elong)) + 
    geom_point(aes(color=factor(biotype))) + 
    facet_grid(rows = ~biotype) +
    labs(x = "GC% of first 1kb", y = "Rate of Elongation", title = "Elongation Rate vs. GC content, split by biotype") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  print(four)
  
  dev.off()
  
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





