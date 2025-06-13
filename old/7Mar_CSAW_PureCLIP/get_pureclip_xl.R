library(dplyr)
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)
library(xlsx)
library(csaw)
library(stringr)

writebed = function(df, fname, bed6 = T){
  
  if(bed6 == T){
    write.table(df[,c("seqnames","start","end","txname","adj_score","strand")], 
                file = fname,
                quote = F,
                row.names = F,
                col.names = F,
                sep = '\t')
  }else{
    write.table(df[,c("seqnames","start","end","score","adj_score","strand","txname")], 
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
  stop("At least 2 arguments must be supplied - xl file and peak regions", call.=FALSE)
}

xl = import(args[1])
peaks = import(args[2])

hits = findOverlaps(xl,peaks)
score_hits = peaks[subjectHits(hits)]$score
xl_tx_names = peaks[subjectHits(hits)]$name

rel_xl = xl[queryHits(hits)]
rel_xl$adj_score = score_hits
rel_xl$txname = xl_tx_names

writebed(as.data.frame(rel_xl), "xl_sites.bed", bed6 = F)

####unique - one CL (top CL score)  per tx...
rel_xl$cluster = csaw::mergeWindows(rel_xl, tol = 1000, ignore.strand = F)$ids

rel_xl_coll = as.data.frame(rel_xl) %>% 
  dplyr::group_by(cluster) %>% 
  arrange(desc(score), .by_group = T) %>%
  dplyr::filter(row_number()==1)

writebed(rel_xl_coll, "xl_sites_collapsed1kb.bed", bed6 = F)




