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
source("~/enhancedClip/functions_2.R")

#annotation_collapsed = read.table("/Users/IEO5559/prog/t3/annotation/rpkm_longest_transcript_whole.tsv", header=T)
annotation_deseq = read.table("/Users/IEO5559/prog/t3/annotation/rpkm_merged_transcripts_w_diffexp.bed", header=F)
annotation = read.table("/Users/IEO5559/prog/t3/annotation/allreads_merged_w_gencode.sorted.gtf", sep = "\t")
h3k4me3 = import("/Users/IEO5559/Desktop/misc/databases/HCT_116_Histone_marks/H3K4me3.bed6.bed")
h3k27ac = import("/Users/IEO5559/Desktop/misc/databases/HCT_116_Histone_marks/H3K27ac.bed6.bed")
Ser5P = import("/Users/IEO5559/Desktop/misc/databases/HCT_116_Histone_marks/Ser5P.bed6.bed")

bams_f = c("/Volumes/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/70_SE_fwd.bam",
           "/Volumes/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/94_SE_fwd.bam",
           "/Volumes/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/95_SE_fwd.bam")

bams_r = c("/Volumes/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/70_SE_rev.bam",
           "/Volumes/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/94_SE_rev.bam",
           "/Volumes/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/95_SE_rev.bam")

rnaseq_f = readGAlignments("/Volumes/data/GN2/gmandana/take2project/t3_reseq/mapped/30/fwd.bam")
rnaseq_r = readGAlignments("/Volumes/data/GN2/gmandana/take2project/t3_reseq/mapped/30/rev.bam")
#peak_regions = import("/Users/IEO5559/enhancedClip/peak_files_for_rpkm_Feb20/20Feb_CSAW_rpkm.bed6.bed")
#peak_table = read.xlsx("/Users/IEO5559/enhancedClip/peak_files_for_rpkm_Feb20/20Feb_CSAW_rpkm.xlsx", sheetIndex = 1)

annotation_ranges = GRanges(seqnames = annotation$V1,
                            ranges = IRanges(start = annotation$V4, end = annotation$V5),
                            strand = annotation$V7)

### use promoters to resize, then overlap with histone marks

c1 = countOverlaps(promoters(annotation_ranges, upstream = 1, downstream = 1), h3k4me3)
c2 = countOverlaps(promoters(annotation_ranges, upstream = 1, downstream = 1), h3k27ac)
c3 = countOverlaps(promoters(annotation_ranges, upstream = 1, downstream = 1), Ser5P)

#don't want exons, want at least one of H3K4me3 and H3K27ac, and Ser5P Pol II
filter = (c1 > 0 | c2 > 0) & annotation$V3 == "transcript" & c3 > 0
annotation_filt = annotation[filter,]

###expression - link to the collapsed tsv
txids_of_tss = unlist(lapply(str_split(lapply(str_split(annotation_filt$V9, ";"), function(x){trimws(x[2])})," "),
                             function(y){y[2]}))
annotation_filt$txid = txids_of_tss

geneid_of_tss = unlist(lapply(str_split(lapply(str_split(annotation_filt$V9, ";"), function(x){trimws(x[1])})," "),
                             function(y){y[2]}))
annotation_filt$geneid = geneid_of_tss

annotation_filt_merged = merge(annotation_filt, annotation_deseq, by.x = 'geneid', by.y = 'V4')

primary = annotation_filt_merged %>% 
  dplyr::select(-c(V9.x, V2.x, V3.x, V6.x, V8.x, V1.y, V2.y, V3.y, V6.y))

colnames(primary) = c('geneid', 'seqid', 'start','end','strand', 'txid','rpkm', 'baseMean', 'log2FC','log10padj')

### expressed RNA only
summary(log(primary$rpkm+0.1))[2]
##first quartile upwards
primary_expr_filt = primary[log(primary$rpkm+0.1) > summary(log(primary$rpkm+0.1))[2],]

primary_expr_filt_ranges = makeGRangesFromDataFrame(primary_expr_filt, 
                                              ignore.strand = F, 
                                              seqnames.field = "seqid", 
                                              start.field = "start", 
                                              end.field = "end", 
                                              strand.field = "strand", 
                                              keep.extra.columns = T)

##how are we collapsing - want to take ratio between upsteram and downstream regions
to_dedup = primary_expr_filt_ranges[duplicated(primary_expr_filt_ranges$geneid)]
to_dedup_up = promoters(to_dedup, upstream = 1000, downstream = 0)
to_dedup_down = promoters(to_dedup, upstream = 0, downstream = 1000)
to_dedup_up_plus = to_dedup_up[strand(to_dedup_up) == "+"]
to_dedup_up_minus = to_dedup_up[strand(to_dedup_up) == "-"]
to_dedup_down_plus = to_dedup_down[strand(to_dedup_down) == "+"]
to_dedup_down_minus = to_dedup_down[strand(to_dedup_down) == "-"]

count_up_plus = countOverlaps(to_dedup_up_plus, rnaseq_f)
count_up_minus = countOverlaps(to_dedup_up_minus, rnaseq_r)
count_down_plus = countOverlaps(to_dedup_down_plus, rnaseq_f)
count_down_minus = countOverlaps(to_dedup_down_minus, rnaseq_r)

to_dedup_plus = to_dedup[strand(to_dedup)=="+"]
to_dedup_minus = to_dedup[strand(to_dedup)=="-"]
to_dedup_plus$ratio = count_down_plus/(count_up_plus+1)
to_dedup_minus$ratio = count_down_minus/(count_up_minus+1)
to_dedup_plus$downcount = count_down_plus
to_dedup_minus$downcount = count_down_minus

to_dedup_new = c(to_dedup_plus,to_dedup_minus)

dedup = as.data.frame(to_dedup_new) %>% 
  group_by(geneid) %>%
  filter(downcount >= mean(downcount)) %>% ##this put here for low expressed "TSS" upstream
  arrange(desc(ratio), by.group = T) %>%
  filter(row_number() == 1)

#rm(to_dedup, to_dedup_down, to_dedup_down_minus, to_dedup_down_plus, to_dedup_plus)
##back to granges
dedup_ranges = makeGRangesFromDataFrame(dedup, 
                                                    ignore.strand = F, 
                                                    seqnames.field = "seqnames", 
                                                    start.field = "start", 
                                                    end.field = "end", 
                                                    strand.field = "strand", 
                                                    keep.extra.columns = T)


unique = c(primary_expr_filt_ranges[!(primary_expr_filt_ranges$geneid %in% dedup_ranges$geneid)],
           dedup_ranges)
unique = unique[!is.na(unique$log10padj)]

#want pval greater than, let's say 0.001, so -log10pval less than 3
five_prime = as.data.frame(reduce(promoters(unique[unique$log10padj < 5], upstream = 0, downstream = 500)))
five_prime$name = "."; five_prime$score = "0"
write.table(five_prime[c(1,2,3,6,7,5)], "all5prime_expressed_unique_500bp_notdiffexpr.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)

######now for the non-bound ones 
unique_resized = promoters(unique, upstream = 0, downstream = 1000)

prom_plus = unique_resized[strand(unique_resized)=="+"]
prom_minus = unique_resized[strand(unique_resized)=="-"]
param = readParam(minq=1, pe="none")
samples = 1:3
Desc = c("INPUT", "IP2", "IP2")
metadata = data.frame(samples, Desc)
desmat <- model.matrix(~factor(metadata$Desc, levels = c("INPUT", "IP2")))
colnames(desmat) = c("INPUT", "IP2")

promcounts_plus <- regionCounts(bams_f, prom_plus, param=param)
promcounts_minus <- regionCounts(bams_r, prom_minus, param=param)

# y_f = asDGEList(promcounts_plus, assay.id=1)
# y_f = estimateDisp(y_f,desmat)
# fit_f = glmQLFit(y_f, desmat, robust=T)
# results_f = glmQLFTest(fit_f, coef="IP2")
# 
# y_r = asDGEList(promcounts_minus, assay.id=1)
# y_r = estimateDisp(y_r,desmat)
# fit_r = glmQLFit(y_r, desmat, robust=T)
# results_r = glmQLFTest(fit_r, coef="IP2")


###compare the CLIP counts with the rpkm value of the transcript, cos rpkm == expression
###so we want tx with nice rpkm but low CLIP counts, right?
###if we divide the rpkm by the CLIP counts, we could get some measure of this
#saveRDS(promcounts_plus, "promcounts_plus.RDS")
#saveRDS(promcounts_minus, "promcounts_minus.RDS")

prom_plus$clipcounts = assay(promcounts_plus)
prom_minus$clipcounts = assay(promcounts_minus)

prom_plus$test = prom_plus$rpkm/(prom_plus$clipcounts + 1)
prom_minus$test = prom_minus$rpkm/(prom_minus$clipcounts + 1)

proms = c(prom_plus,prom_minus)

####what thresholds should I use to define non-bound peaks?
summary(as.vector(proms$test))
summary(as.vector(proms$clipcounts))
test = proms[proms$test > 0.010469]

final_filter = as.vector(proms$test > median(proms$test) & proms$clipcounts < median(proms$clipcounts))

five_prime_bg2 = as.data.frame(reduce(promoters(proms[final_filter], upstream = 0, downstream = 500)))
five_prime_bg2$name = "."; five_prime_bg2$score = "0"
write.table(five_prime_bg2[c(1,2,3,6,7,5)], "all5prime_expressed_unique_500bp_notdiffexpr_notbound.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)



