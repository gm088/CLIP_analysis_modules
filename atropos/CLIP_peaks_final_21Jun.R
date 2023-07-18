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
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")

bams_f = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_fwd_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/94_SE_fwd_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/95_SE_fwd_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105289/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105291/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105293/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105295/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105290/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105292/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105294/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105296/fwd.bam")

bams_r = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_rev_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/94_SE_rev_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/95_SE_rev_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105289/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105291/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105293/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105295/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105290/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105292/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105294/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/4sU_Danilo_HCT116/BAM/SID105296/rev.bam")

samples = 1:11
Desc = c("INPUT", "IP2", "IP2", "CTL", "CTL", "CTL", "CTL", "AUX", "AUX", "AUX", "AUX")
metadata = data.frame(samples, Desc)
desmat <- model.matrix(~factor(metadata$Desc, levels = c("INPUT", "IP2", "CTL", "AUX")))
colnames(desmat) = c("INPUT", "IP2", "CTL", "AUX")

param = readParam(minq=10, pe="both")

#this is viviana's
new_siW = customimportbed6('new_siW/all_siW_new.bed')
peaks_file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_20Jun/clipper_all_rpkm.bed"
peaks_file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_20Jun/clipper_new_rpkm_noucscrepeats_nosmallRNA.bed"
peaks = customimportbed6(peaks_file)
colnames(mcols(peaks)) = c("txname",
                          "adj_score",
                          "log10pval",
                          "log2FC",
                          "rpkm")

####pivot
clusters = as.data.frame(peaks) %>% 
  group_by(txname) %>% 
  mutate(cluster_start = min(start)) %>% 
  mutate(cluster_end = max(end)) %>% 
  mutate(coords = paste(seqnames, paste(start, end, sep = "-"), sep = ":")) %>% 
  select(-c("start", "end", "width")) %>% 
  summarise_all(paste, collapse = ",")

#for later
mean_scores = as.data.frame(peaks) %>% 
  group_by(txname) %>% 
  summarise(cluster_mean_score = mean(adj_score))

clusters$seqnames = unlist(lapply(clusters$seqnames, 
                                                  function(x){str_split(x, ",")[[1]][[1]]}))
clusters$strand = unlist(lapply(clusters$strand, 
                                  function(x){str_split(x, ",")[[1]][[1]]}))
clusters$cluster_start = as.numeric(unlist(lapply(clusters$cluster_start, 
                                                  function(x){str_split(x, ",")[[1]][[1]]})))
clusters$cluster_end = as.numeric(unlist(lapply(clusters$cluster_end, 
                                                function(x){str_split(x, ",")[[1]][[1]]})))

clusters = clusters[,c(2,8,9,1,3,4,5,6,7,10)]

clusters_ranges = makeGRangesFromDataFrame(clusters,
                                          start.field = "cluster_start",
                                          end.field = "cluster_end",
                                          strand.field = "strand",
                                          keep.extra.columns = T)

clus_olap = findOverlaps(clusters_ranges, new_siW, maxgap = 500, ignore.strand = F)

##append the overlap info
clusters[queryHits(clus_olap)[!duplicated(queryHits(clus_olap))],]

toappend = as.data.frame(mcols(new_siW[subjectHits(clus_olap)[!duplicated(queryHits(clus_olap))]])) %>% 
  mutate(summ = paste(V4, V7, V8, V9, sep = ";"))

clusters$overlap = NA
clusters[queryHits(clus_olap)[!duplicated(queryHits(clus_olap))],]$overlap = toappend$summ

clusters = as_tibble(merge(clusters, mean_scores))

###count at the cluster level

regions_plus = resize(clusters_ranges[strand(clusters_ranges)=="+"], 
                      width = 1000, fix = "center")
regions_minus = resize(clusters_ranges[strand(clusters_ranges)=="-"],
                       width = 1000, fix = "center")

#count for and rev separately - count ONLY RNA SEQ BAMS, check param
regcounts_plus <- regionCounts(bams_f[4:11], regions_plus, param=param)
regcounts_minus <- regionCounts(bams_r[4:11], regions_minus, param=param)

regions_plus$ctlcounts = rowMeans(assay(regcounts_plus)[,c(1:4)])
regions_plus$auxcounts = rowMeans(assay(regcounts_plus)[,c(5:8)])

regions_minus$ctlcounts = rowMeans(assay(regcounts_minus)[,c(1:4)])
regions_minus$auxcounts = rowMeans(assay(regcounts_minus)[,c(5:8)])

regions_all = c(regions_plus, regions_minus)

clusters = as_tibble(merge(clusters, as.data.frame(regions_all)[,c(6,12,13)]))

clusters_final = clusters[order(clusters$cluster_mean_score, decreasing = T),]

saveRDS(clusters_final, file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_20Jun/clusters_final.RDS")

#### consolidated peaks, individual
peaks_table = as.data.frame(peaks)

peaks_olap = findOverlaps(peaks, new_siW, maxgap = 500, ignore.strand = F)

##append the overlap info
peaks_table[queryHits(peaks_olap)[!duplicated(queryHits(peaks_olap))],]

toappend = as.data.frame(mcols(new_siW[subjectHits(peaks_olap)[!duplicated(queryHits(peaks_olap))]])) %>% 
  mutate(summ = paste(V4, V7, V8, V9, sep = ";"))

peaks_table$overlap = NA
peaks_table[queryHits(peaks_olap)[!duplicated(queryHits(peaks_olap))],]$overlap = toappend$summ

peaks_table = as_tibble(merge(peaks_table, mean_scores))

saveRDS(peaks_table, file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_20Jun/idr_peaks_5Jul_nosmallRNA.RDS")


#####rep1 and rep2

stem = "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/clipper/using_siW_w_coding/nomismatch/outputs/21Jun/"

rep1info = read.table(paste0(stem, "rep1_sig.bed.full"), sep = "\t", header = F)
rep2info = read.table(paste0(stem, "rep2_sig.bed.full"), sep = "\t", header = F)

rep1ranges = customimportbed6(paste0(stem, "rep1_sig.bed"))
rep2ranges = customimportbed6(paste0(stem, "rep2_sig.bed"))

###overlap with siW - rep1

rep1_olap = findOverlaps(rep1ranges, new_siW, maxgap = 500, ignore.strand = F)

##append the overlap info
rep1info[queryHits(rep1_olap)[!duplicated(queryHits(rep1_olap))],]

toappend = as.data.frame(mcols(new_siW[subjectHits(rep1_olap)[!duplicated(queryHits(rep1_olap))]])) %>% 
  mutate(summ = paste(V4, V7, V8, V9, sep = ";"))

rep1info$overlap = NA
rep1info[queryHits(rep1_olap)[!duplicated(queryHits(rep1_olap))],]$overlap = toappend$summ

rep1info_final = as.data.frame(rep1info) %>% 
  select(-c("V4", "V7", "V8", "V9", "V10"))

colnames(rep1info_final) = c("seqnames",
                             "start",
                             "end",
                             "IP_counts",
                             "SMInput_counts",
                             "log10pval",
                             "log2FC",
                             "overlap_with_ZC3H4_sensitive")

rep1info_final$strand = as.vector(strand(rep1ranges))
#add peakID
rep1info_final$peak_id = paste0("id", 1:nrow(rep1info_final))
saveRDS(rep1info_final, file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_20Jun/rep1info_final.RDS")

#### rep 2

rep2_olap = findOverlaps(rep2ranges, new_siW, maxgap = 500, ignore.strand = F)

##append the overlap info
rep2info[queryHits(rep2_olap)[!duplicated(queryHits(rep2_olap))],]

toappend = as.data.frame(mcols(new_siW[subjectHits(rep2_olap)[!duplicated(queryHits(rep2_olap))]])) %>% 
  mutate(summ = paste(V4, V7, V8, V9, sep = ";"))

rep2info$overlap = NA
rep2info[queryHits(rep2_olap)[!duplicated(queryHits(rep2_olap))],]$overlap = toappend$summ

rep2info_final = as.data.frame(rep2info) %>% 
  select(-c("V4", "V7", "V8", "V9", "V10"))

colnames(rep2info_final) = c("seqnames",
                             "start",
                             "end",
                             "IP_counts",
                             "SMInput_counts",
                             "log10pval",
                             "log2FC",
                             "overlap_with_ZC3H4_sensitive")

rep2info_final$strand = as.vector(strand(rep2ranges))

rep2info_final$peak_id = paste0("id", 1:nrow(rep2info_final))

#export RDS for excel writng
saveRDS(rep2info_final, file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_20Jun/rep2info_final.RDS")

##homogenise the peaks

newtop = peaks_new[!duplicated(peaks_new$name)]

newtop$quintile = cut(newtop$score, 
                   breaks = rev(quantile(newtop$score, probs = seq(0,1,0.2))), 
                   labels = c(1,2,3,4,5), 
                   include.lowest = T)

newtop$quintile = factor(newtop$quintile, levels = c(5,4,3,2,1))

unlist(lapply(5:1, function(x){length(subsetByOverlaps(newtop[newtop$quintile == x], 
                                                  siW))/length(newtop[newtop$quintile == x])}))


##### 20 Jun - the initial peakcalling, but with histones etc. and not just top 1000

proms = promoters(new_peaks, upstream = 0, downstream = 1000)
regions2 = proms[!duplicated(csaw::mergeWindows(proms, tol = 0)$ids)]
regions2$id = paste("id", 1:length(regions2), sep = "")

regions_plus2 = regions2[strand(regions2)=="+"]
regions_minus2 = regions2[strand(regions2)=="-"]

#count for and rev separately - count ONLY RNA SEQ BAMS, check param
regcounts_plus <- regionCounts(bams_f[4:9], regions_plus2, param=param)
regcounts_minus <- regionCounts(bams_r[4:9], regions_minus2, param=param)

regions_plus2$ctlcounts = log2(rowMeans(assay(regcounts_plus)[,c(1:3)]))
regions_plus2$auxcounts = log2(rowMeans(assay(regcounts_plus)[,c(4:6)]))

regions_minus2$ctlcounts = log2(rowMeans(assay(regcounts_minus)[,c(1:3)]))
regions_minus2$auxcounts = log2(rowMeans(assay(regcounts_minus)[,c(4:6)]))

peaks_new2 = c(regions_plus2, regions_minus2)
peaks_new2 = peaks_new2[order(peaks_new2$V5, decreasing = T)]

### quintiles of peak scores, with overlap with siW, with violin plots
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

peaks_new2$quintile = cut(peaks_new2$V5, 
                          breaks = rev(quantile(peaks_new2$V5, probs = seq(0,1,0.2))), 
                          labels = c(1,2,3,4,5), 
                          include.lowest = T)

peaks_new2$quintile = factor(peaks_new2$quintile, levels = c(5,4,3,2,1))

ol = unlist(lapply(5:1, function(x){length(subsetByOverlaps(peaks_new2[peaks_new2$quintile == x], new_siW, maxgap = 500))}))

as.data.frame(peaks_new2) %>% 
  pivot_longer(., cols = ends_with("counts"), names_repair = "unique") %>%
  mutate_at("name", as.factor) %>%
  mutate(name = fct_relevel(name, "ctlcounts", "auxcounts")) %>%
  ggplot(aes(quintile, value)) + 
  geom_violin(aes(fill = name), width=0.6, alpha = 0.4) + 
  geom_boxplot(aes(fill = name), width=0.1, color="grey", alpha=0.8, position=position_dodge(0.6)) + 
  theme_minimal() + 
  scale_fill_manual(values=cbbPalette) + 
  xlab("")

colnames(mcols(peaks_new2)) = c("txid",
                                "adj_score",
                                "log2FC",
                                "log10pval",
                                "rpkm",
                                "id",
                                "ctlcounts",
                                "auxcounts",
                                "quintile")




write.table(as.data.frame(regions2)[,c(1,2,3,6,7,5,8,9,10,12)],
            file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_20Jun/peaks_1kbext.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

write.table(as.data.frame(new_peaks)[,c(1,2,3,6,7,5,8,9, 10,11)],
            file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_20Jun/peaks_w_quintile.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
















