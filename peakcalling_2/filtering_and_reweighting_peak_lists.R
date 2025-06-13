library(dplyr)
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)
library(xlsx)
library(rtracklayer)

################

## this script was made to summarise the final filering and rpkm weighting steps
## especially because I needed to process the relaxed IDR files

##look to rpkm_weighting for original script

###############

customimportbed6 = function(fname){
  
  regions = makeGRangesFromDataFrame(read.table(fname, header = F),
                                     keep.extra.columns = T, ignore.strand = F,
                                     seqnames.field = "V1",start.field = "V2",
                                     end.field = "V3",strand.field = "V6")
  
  return(regions)
}

getann = function(peaksfile, annot, handleNA = T){
  
  hits = findOverlaps(peaksfile, annot)
  unann = subsetByOverlaps(peaksfile, annot, invert = T)
  int1 = peaksfile[queryHits(hits)]
  #append the annot info the overlapping peaks
  int1$txname = annot[subjectHits(hits)]$name
  int1$rpkm = annot[subjectHits(hits)]$score
  #add the unann peaks
  if(handleNA == T & length(unann) > 0){
    unann$rpkm = 0
    final = c(int1, unann)
    final$adj_score = final$score/(final$rpkm + 0.1)
    return(return(final[order(seqnames(final))]))
  }else{
    final = int1
    final$adj_score = final$score/(final$rpkm + 0.1)
    return(final[order(seqnames(final))])
  }
}

adjust = function(df){
  
  #this is to account for one unique peak being on two transcripts
  df1 = as.data.frame(df) %>% 
    mutate(coor = paste(as.character(seqnames), as.character(start), as.character(end), sep = ".")) %>% 
    dplyr::group_by(coor) %>% 
    arrange(adj_score) %>% 
    dplyr::filter(row_number()==1) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-coor)
  return(df1[order(df1$adj_score, decreasing = T),])
}

process_file = function(filepath, minwidth, annotation_file = "/Users/IEO5559/enhancedClip/10may_ttome/rpkm_merged_transcripts.bed",
                        rmrpts = F, write = F, outdir, fname, ext_ann = T, handleNAs = F){
  
  #preproc
  peaksfile = customimportbed6(filepath)
  annotation = customimportbed6(annotation_file)
  if(ext_ann == T){
    annotation = resize(annotation, fix = "center", width = width(annotation)+2000)
  }
  colnames(mcols(annotation)) = c("name", "score")
  colnames(mcols(peaksfile)) = c("score", "name")
  peaksfile$score = as.numeric(peaksfile$score)
  #browser()
  head(peaksfile)
  head(annotation)
  
  ###tracking
  peaksfile$init_id = 1:length(peaksfile)
  
  #smallRNAtothrow
  smallRNA = customimportbed6("/Users/IEO5559/Desktop/misc/annotations/small_RNA_tRNA_RPL_RPPH_bed6.bed")
  peaks = subsetByOverlaps(peaksfile, smallRNA, maxgap = 100, invert = T)
  
  ###this is for tRNA
  if(rmrpts == T){
    peaks = subsetByOverlaps(peaks, 
                             customimportbed6("/Users/IEO5559/Desktop/misc/databases/ucsc/justRNA.bed"), 
                             invert = T)
  }
  
  ##too short
  peaks = peaks[width(peaks) > minwidth,]
  
  #rpkm weighting
  clipper_rpkm = adjust(getann(peaks, annotation, handleNA = handleNAs))
  
  if(write == T){
    
    write.table(clipper_rpkm[,c("seqnames", "start", "end", "txname", "adj_score", "strand", "score", "name", "rpkm")],
                file = paste0(outdir, fname),
                quote = F,
                col.names = F,
                row.names = F,
                sep = "\t")
    
  }
  
  return(clipper_rpkm)
}

####### 6 oct - ZC3H4
dir = "/Users/IEO5559/enhancedClip/ip7_reseq/"

process_file(paste0(dir, "reproducible_peaks_sorted.bed"), 
             rmrpts = T, write = T, outdir = dir, 
             fname = "reproducible_peaks_reseq_rpkm.bed")


###### DRB 5 min, siW3_refined_ttome
dir = "/Users/IEO5559/enhancedClip/DRB_5min_siW3_refined/"
peaks_adjusted = process_file(paste0(dir, "reproducible_peaks_sorted.bed"), 
                              annotation_file = "~/prog/DRB_TTseq/siW_ttome/t3_quant/siW3_refined_merged_transcripts_t3.bed",
                              rmrpts = T, write = T, outdir = dir, 
                              fname = "DRB5_idr_peaks_rpkm_t3.bed", handleNAs = F, 
                              ext_ann = T, minwidth = 30)













