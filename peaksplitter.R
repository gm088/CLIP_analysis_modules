library(csaw)
library(rtracklayer)
library(edgeR)
library(BiocParallel)
library(GenomicRanges)
library(GenomicAlignments)
library(slider)
source("utility_functions.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")

######  Thurs 15th Dec
###load bamfiles for making coverage
bams_f = c("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/IP7_DCIP_SE_fwd.bam")

bams_r = c("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/IP7_DCIP_SE_rev.bam")

### Use with_clean_bams_12Dec_2022_no_tRNA_final_win150bp_spacing50bp.bed to define peak regions

#regions = import("with_clean_bams_12Dec_2022_no_tRNA_final_win150bp_spacing50bp.bed")
#bigwin = import("with_clean_bams_6Dec_2022_no_tRNA_final_win500bp.bed")
#reg_clean = subsetByOverlaps(regions, bigwin, ignore.strand = F)
#names(reg_clean) = 1:length(reg_clean)

regions = import("with_clean_bams_6Dec_2022_no_tRNA_final_win500bp.bed")
names(regions) = paste("id", 1:length(regions), sep = "")

r_plus = regions[strand(regions)=="+"]
r_minus = regions[strand(regions)=="-"]

#careful, this function requires your bedfile to have names i.e names(regions) must not be NULL

plus_prof = coverageFromFile(txdb = NULL, file = bams_f, gr = r_plus, change_chr_names = F, strand = "+")
minus_prof = coverageFromFile(txdb = NULL, file = bams_r, gr = r_minus, change_chr_names = F, strand = "-")

###signal processing

peaksplitter = function(vector, winsize = 10, thresh = 0.07){
  
  test = unlist(slide(vector, mean, .before = winsize, .after = winsize))
  ddxtest = test[2:length(test)] - test[1:(length(test)-1)]
  ddx2test = ddxtest[2:length(ddxtest)] - ddxtest[1:(length(ddxtest)-1)]
  
  ddxtest_smooth = unlist(slide(ddxtest, mean, .before = winsize, .after = winsize))
  ddx2test_smooth = unlist(slide(ddx2test, mean, .before = winsize, .after = winsize))
  max = slide_max(ddxtest, before = winsize, after = winsize)
  min = slide_min(ddxtest, before = winsize, after = winsize)
  
  indices = max > thresh & min < -thresh
  peaks = split(which(indices), cumsum(seq_along(which(indices)) %in% (which(diff(which(indices))>1)+1)))
  
  #additional - just make sure all are maxima
  peaks_filt = peaks[unlist(lapply(peaks, function(x){sum(ddx2test[x]) < 0}))]
  
  #this filtering step cos there were some subpeaks just listed as "NULL"
  peaks_filt = lapply(peaks_filt, function(x){x[lapply(x, length) > 0]})
  return(peaks_filt)
}

break_them_up = function(interval, indices, strand){
  
  chr = as.character(seqnames(interval))
  
  start = start(interval)
  newstarts = unlist(lapply(indices, function(x){start + x[1]}))
  newends = unlist(lapply(indices, function(x){start + x[length(x)]}))
  all_split_up = GRanges(seqnames = rep(chr, length(indices)),
                         ranges = IRanges(start = newstarts, end = newends),
                         strand = rep(strand, length(indices)))
  
  return(all_split_up)
}

adjust = function(granges, splitpeaks, strand){
  
  stopifnot(length(granges) == length(splitpeaks))
  new = lapply(1:length(granges), function(x){break_them_up(granges[x], splitpeaks[[x]], strand)})
  return(new)
  
}

check_plus = lapply(plus_prof, function(x){peaksplitter(x)})
ready2split_plus = lapply(check_plus, function(x){x[lapply(x, length)>0]})
split_plus = adjust(r_plus, ready2split_plus, "+")

check_minus = lapply(minus_prof, function(x){peaksplitter(x)})
ready2split_minus = lapply(check_minus, function(x){x[lapply(x, length)>0]})
split_minus = adjust(r_minus, ready2split_minus, "-")

testoutput = GRangesList(c(split_plus, split_minus))

checkme = export(unlist(testoutput), "split_peaks_IP7DCIP_19Dec22.bed")



















