library(csaw)
library(rtracklayer)
library(edgeR)
library(BiocParallel)
library(GenomicRanges)
library(GenomicAlignments)
library(slider)
library(dplyr)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(fitdistrplus)
source("~/enhancedClip/functions_2.R")

annotation_collapsed = read.table("/Users/IEO5559/prog/t3/annotation/rpkm_longest_transcript_whole.tsv")
annotation = read.table("/Users/IEO5559/prog/t3/annotation/allreads_merged_w_gencode.sorted.gtf", sep = "\t")
h3k4me3 = import("/Users/IEO5559/Desktop/misc/databases/HCT_116_Histone_marks/H3K4me3.bed6.bed")
h3k27ac = import("/Users/IEO5559/Desktop/misc/databases/HCT_116_Histone_marks/H3K27ac.bed6.bed")

peak_regions = import("/Users/IEO5559/enhancedClip/peak_files_for_rpkm_Feb20/20Feb_CSAW_rpkm.bed6.bed")
peak_table = read.xlsx("/Users/IEO5559/enhancedClip/peak_files_for_rpkm_Feb20/20Feb_CSAW_rpkm.xlsx", sheetIndex = 1)

annotation_ranges = GRanges(seqnames = annotation$V1,
                            ranges = IRanges(start = annotation$V4, end = annotation$V5),
                            strand = annotation$V7)
###use H3K4me3 and H3K27ac to determine start sites
### use promoters to resize, then overlap with histone marks

c1 = countOverlaps(promoters(annotation_ranges, upstream = 1000, downstream = 1000), h3k4me3)
c2 = countOverlaps(promoters(annotation_ranges, upstream = 1000, downstream = 1000), h3k27ac)

#don't want exons
filter = (c1 > 0 & c2 > 0) & annotation$V3 == "transcript"
annotation$istss = filter
annotation_filt = annotation[filter,]

###expression - link to the collapsed tsv
txids_of_tss = unlist(lapply(str_split(lapply(str_split(annotation_filt$V9, ";"), function(x){trimws(x[2])})," "),
                      function(y){y[2]}))

tx_main_tss = annotation_collapsed[annotation_collapsed$V7 %in% txids_of_tss,]

###not a peak region
tx_main_tss_ranges = makeGRangesFromDataFrame(tx_main_tss, 
                         ignore.strand = F, 
                         seqnames.field = "V1", 
                         start.field = "V2", 
                         end.field = "V3", 
                         strand.field = "V4", 
                         keep.extra.columns = T)

tx_main_tss_nopeak = subsetByOverlaps(tx_main_tss_ranges, peak_regions, invert = T)
tx_main_tss_nopeak$V5 = as.numeric(tx_main_tss_nopeak$V5)

###expression - we want it to be of roughly same rpkm
summary(peak_table$rpkm)
summary(tx_main_tss_nopeak$V5)

par(mfrow = c(1,2))
hist(log10(peak_table$rpkm), breaks = seq(-4,4,by = 0.1))
hist(log10(tx_main_tss_nopeak$V5), breaks = seq(-4,4,by = 0.1))

###based on these histograms let's say take tx with log10(rpkm) between -1 and 1
min_log10_rpkm = -1
max_log10_rpkm = 1

tx_main_tss_nopeak_exprcontrol = tx_main_tss_nopeak[log10(tx_main_tss_nopeak$V5) > min_log10_rpkm & 
                                                      log10(tx_main_tss_nopeak$V5) < max_log10_rpkm]


####### this is the "final" list of tx which have TSS WITHOUT ENRICHMENT
tx_main_tss_nopeak_exprcontrol

### take promoter regions
promoters(tx_main_tss_nopeak_exprcontrol, upstream = 400, downstream = 100)
### randomly sample 1000
sample(promoters(tx_main_tss_nopeak_exprcontrol, upstream = 400, downstream = 100), 1000)

### write output
export(promoters(tx_main_tss_nopeak_exprcontrol, upstream = 400, downstream = 100),
       "~/enhancedClip/tx_5prime_nopeak_exprcontrol.bed")

export(sample(promoters(tx_main_tss_nopeak_exprcontrol, upstream = 400, downstream = 100), 1000), 
       "~/enhancedClip/tx_5prime_nopeak_exprcontrol_randomsample.bed")









