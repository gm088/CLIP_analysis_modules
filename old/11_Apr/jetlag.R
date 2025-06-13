library(csaw)
library(rtracklayer)
library(edgeR)
library(BiocParallel)
library(corrplot)
library(GenomicAlignments)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(forcats)

bams_f = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_fwd_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/94_SE_fwd_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/95_SE_fwd_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/07/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/15/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/23/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/08/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/16/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/24/fwd.bam")

bams_r = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_rev_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/94_SE_rev_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/95_SE_rev_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/07/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/15/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/23/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/08/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/16/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/24/rev.bam")

samples = 1:9
Desc = c("INPUT", "IP2", "IP2", "CTL", "CTL", "CTL", "AUX", "AUX", "AUX")
metadata = data.frame(samples, Desc)
desmat <- model.matrix(~factor(metadata$Desc, levels = c("INPUT", "IP2", "CTL", "AUX")))
colnames(desmat) = c("INPUT", "IP2", "CTL", "AUX")

param = readParam(minq=10, pe="both")
caller = "clipper"
siW = import('/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mysiW/mysiW_bed6.bed')

if(caller == "csaw"){
  regions = import("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/csaw/14Feb_locfilt_gradfilt_nohistone.bed")
  #because csaw regions width is crazy
  #regions = resize(regions, width = 1000, fix = "start")
}else if(caller == "pchu"){
  regions = import("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/peakachu/peakachu_2reps_nomismatch_mad0_fc3/peaks_peaks.bed")
}else if(caller == "clipper"){
  peaks = import("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/22Feb_rpkm_norm/20Feb_CLIPper_rpkm.bed6.bed")
  top = peaks[1:1000]
}

regions = promoters(top, upstream = 0, downstream = 1000)
regions$id = paste("id", 1:length(regions), sep = "")

regions_plus = regions[strand(regions)=="+"]
regions_minus = regions[strand(regions)=="-"]

#count for and rev separately - count ONLY RNA SEQ BAMS, check param
regcounts_plus <- regionCounts(bams_f[4:9], regions_plus, param=param)
regcounts_minus <- regionCounts(bams_r[4:9], regions_minus, param=param)


########test

#divide counts by width and make distribution
mcols(regions_plus) = cbind(mcols(regions_plus), "countperkb" = (rowMeans(assay(regcounts_plus))/width(regions_plus))*1000)
mcols(regions_minus) = cbind(mcols(regions_minus), "countperkb" = (rowMeans(assay(regcounts_minus))/width(regions_minus))*1000)

regions_new = c(regions_plus, regions_minus)
regions_new_sorted = regions_new[order(regions_new$countperkb)]

###write files
write.table(as.data.frame(regions_new_sorted)[,c(1,2,3,7,8,5)], 
            "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/csaw/14Feb_locfilt_gradfilt_nohistone_sortedbyRNA.bed",
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F)


########## nascent RNA counts

counts = rbind(assay(regcounts_plus), assay(regcounts_minus))
libsizes = regcounts_plus$totals + regcounts_minus$totals
norm_counts = sweep(counts, 2, libsizes, FUN = "/") * 10000000

log_counts = log((norm_counts), base = 2)

wt_rna_abun = rowMeans(log_counts[,c(4,5,6)])
summary(wt_rna_abun)

cut(wt_rna_abun, breaks=c(quantile(wt_rna_abun, probs = seq(0, 1, by = 0.1))), labels=1:10, include.lowest=TRUE)

regions$quantile = cut(wt_rna_abun, breaks=c(quantile(wt_rna_abun, probs = seq(0, 1, by = 0.1))), labels=1:10, include.lowest=TRUE)

par(mfrow=c(2,5))

for(i in 1:10){
  boxplot(norm_counts[regions$quantile==i,], outline = F, 
          col = c('red','blue', 'green','brown', 'brown', 'brown','orange', 'orange', 'orange'),
          main = paste("decile no.", i, sep = ""), xaxt="n")
}

par(mfrow=c(2,5))

for(i in 1:10){
  x = barplot(table(regions[regions$quantile==i]$annot), xaxt="n")
  labs <- names(table(regions[regions$quantile==i]$annot))
  text(cex=0.9, x=x-.25, y=-7, labs, xpd=TRUE, srt=45)
}


######## 10th April 2023 #########

regions_plus$ctlcounts = log(rowMeans(assay(regcounts_plus)[,c(1:3)]))
regions_plus$auxcounts = log(rowMeans(assay(regcounts_plus)[,c(4:6)]))

regions_minus$ctlcounts = log(rowMeans(assay(regcounts_minus)[,c(1:3)]))
regions_minus$auxcounts = log(rowMeans(assay(regcounts_minus)[,c(4:6)]))

peaks_new = c(regions_plus, regions_minus)
peaks_new = peaks_new[order(peaks_new$score, decreasing = T)]

### quintiles of peak scores, with overlap with siW, with violin plots
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

peaks_new$quintile = cut(peaks_new$score, 
                     breaks = rev(quantile(peaks_new$score, probs = seq(0,1,0.2))), 
                     labels = c(1,2,3,4,5), 
                     include.lowest = T)

peaks_new$quintile = factor(peaks_new$quintile, levels = c(5,4,3,2,1))

ol = unlist(lapply(5:1, function(x){length(subsetByOverlaps(peaks_new[peaks_new$quintile == x], siW, maxgap = 500))}))

as.data.frame(peaks_new)[,c(8,9,10,11)] %>% 
  pivot_longer(., cols = ends_with("counts"), names_repair = "unique") %>%
  mutate_at("name", as.factor) %>%
  mutate(name = fct_relevel(name, "ctlcounts", "auxcounts")) %>%
  ggplot(aes(quintile, value)) + 
  geom_violin(aes(fill = name), width=0.6, alpha = 0.4) + 
  geom_boxplot(aes(fill = name), width=0.1, color="grey", alpha=0.8, position=position_dodge(0.6)) + 
  theme_minimal() + 
  scale_fill_manual(values=cbbPalette) + 
  xlab("")

#####ones that overlap with siW

peaks_siW = subsetByOverlaps(peaks_new, siW, maxgap = 500)
sample_size = as.data.frame(peaks_siW) %>% group_by(quintile) %>% summarize(num=n())

as.data.frame(peaks_siW)[,c(8,9,10,11)] %>% 
  left_join(sample_size) %>% 
  pivot_longer(., cols = ends_with("counts"), names_repair = "unique") %>%
  mutate_at("name", as.factor) %>%
  mutate(name = fct_relevel(name, "ctlcounts", "auxcounts")) %>%
  mutate(myaxis = factor(paste0(quintile, ", ", "n=", num), levels = unique(paste0(quintile, ", ", "n=", num)))) %>%
  ggplot(aes(myaxis, value)) + 
  geom_violin(aes(fill = name), width=0.6, alpha = 0.4) + 
  geom_boxplot(aes(fill = name), width=0.1, color="grey", alpha=0.8, position=position_dodge(0.6)) + 
  theme_minimal() + 
  scale_fill_manual(values=cbbPalette) + 
  xlab("")

#####ones that dont overlap with siW

peaks_no_siW = subsetByOverlaps(peaks_new, siW, maxgap = 500, invert = T)
sample_size = as.data.frame(peaks_no_siW) %>% group_by(quintile) %>% summarize(num=n())

as.data.frame(peaks_no_siW)[,c(8,9,10,11)] %>% 
  left_join(sample_size) %>% 
  pivot_longer(., cols = ends_with("counts"), names_repair = "unique") %>%
  mutate_at("name", as.factor) %>%
  mutate(name = fct_relevel(name, "ctlcounts", "auxcounts")) %>%
  mutate(myaxis = factor(paste0(quintile, ", ", "n=", num), levels = unique(paste0(quintile, ", ", "n=", num)))) %>%
  ggplot(aes(myaxis, value)) + 
  geom_violin(aes(fill = name), width=0.6, alpha = 0.4) + 
  geom_boxplot(aes(fill = name), width=0.1, color="grey", alpha=0.8, position=position_dodge(0.6)) + 
  theme_minimal() + 
  scale_fill_manual(values=cbbPalette) + 
  xlab("")

######characterising peaks that don't overlap with siW

###my crutches

ser5 = import("/hpcnfs/data/GN2/gmandana/annotation/Ser5P.bed6.bed")
k27ac = import("/hpcnfs/data/GN2/gmandana/annotation/H3K27ac.bed6.bed")
enh = import("/hpcnfs/data/GN2/gmandana/annotation/hct116enh.bed")

ann_deseq2 = read.table("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/stringtie/rpkm_merged_transcripts_w_diffexp.bed", header = F)
colnames(ann_deseq2) = c("chr", "start", "stop", "name", "RPKM", "strand", "baseMean", "log2FC", "logpadj")
rownames(ann_deseq2) = ann_deseq2$name

##these need to be defined
top$quintile = cut(top$score, 
                   breaks = rev(quantile(top$score, probs = seq(0,1,0.2))), 
                   labels = c(1,2,3,4,5), 
                   include.lowest = T)

top$quintile = factor(top$quintile, levels = c(5,4,3,2,1))

not_siW_binding = subsetByOverlaps(top, siW, maxgap = 500, invert = T)
not_siW_binding$h3k27ac = countOverlaps(not_siW_binding, k27ac, maxgap = 500) > 0
not_siW_binding$ser5 = countOverlaps(not_siW_binding, ser5, maxgap = 500) > 0

notk27ac = subsetByOverlaps(not_siW_binding, k27ac, maxgap = 500, invert = T)

summary = as.data.frame(not_siW_binding) %>%
  group_by(quintile) %>%
  summarise(sum(as.numeric(h3k27ac)), sum(as.numeric(ser5)))

totals = as.data.frame(not_siW_binding) %>%
  group_by(quintile) %>%
  tally()
