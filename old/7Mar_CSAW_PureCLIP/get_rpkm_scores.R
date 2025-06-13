library(dplyr)
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)
library(xlsx)
library(stringr)

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

writebed = function(df, fname, bed6 = T){
  
  if(bed6 == T){
    write.table(df[,c("seqnames","start","end","txname","adj_score","strand")], 
                file = fname,
                quote = F,
                row.names = F,
                col.names = F,
                sep = '\t')
  }else{
    write.table(df[,c("seqnames","start","end","txname","adj_score","strand","rpkm","width")], 
                file = fname,
                quote = F,
                row.names = F,
                col.names = F,
                sep = '\t')
  }
  #write.xlsx(df, file = paste(outdir, fname, ".xlsx", sep = ""))
  
}

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

outfile = paste(str_split(args[1], "\\.")[[1]][1], "_rpkm.bed", sep = "")

###load the file with rpkms
annotation = import("/Users/IEO5559/prog/t3/annotation/rpkm_merged_transcripts.bed")
peaks = import(args[1])
#return dataframes
peaks_rpkm = adjust(getann(peaks, annotation, handleNA = T))
#write
writebed(peaks_rpkm, outfile)




