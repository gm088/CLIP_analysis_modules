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

bams_f = c("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/07/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/15/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/23/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/08/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/16/fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/24/fwd.bam")

bams_r = c("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/07/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/15/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/23/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/08/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/16/rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/24/rev.bam")

param = readParam(minq=10, pe="both")

#this is mine - refined
siW = readRDS("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/T3_ghostbusters/ZC3H4_supp_tx_PROCAP_FANTOM.RDS")

names(mcols(siW)) = c("name", "FDR", "class", "WT_rpkm", "adj_FDR", "class2")
siW$id = paste0("id", 1:length(siW))

peaks_file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_26Oct/reproducible_peaks_reseq_rpkm.bed"
peaks = customimportbed6(peaks_file)
colnames(mcols(peaks)) = c("txname",
                           "adj_score",
                           "log10pval",
                           "log2FC",
                           "rpkm")

#first one is tower; how did it get there?
peaks[1]$adj_score = 1; peaks[1]$log10pval = 0.1; peaks[1]$log2FC = 3; peaks[1]$rpkm = 0

macspeaks = customimportbed6("/hpcnfs/data/GN2/gmandana/take2project/PROSEQ/macs/round3/macspeaks.bed")
proseqtss = distanceToNearest(peaks, resize(macspeaks, width = 1, fix = "start"))

peaks$proseqtss = (as.data.frame(macspeaks[subjectHits(proseqtss)]) %>%
                     mutate(coords = paste(seqnames, paste(start, end, sep = "-"), sep = ":")) %>% 
                     dplyr::select("coords"))$coords

peaks$dsisttoproseq = mcols(proseqtss)$distance

####pivot
clusters = as.data.frame(peaks) %>% 
  group_by(txname) %>% 
  mutate(cluster_start = min(start)) %>% 
  mutate(cluster_end = max(end)) %>% 
  mutate(coords = paste(seqnames, paste(start, end, sep = "-"), sep = ":")) %>% 
  dplyr::select(-c("start", "end", "width")) %>% 
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

clus_olap = findOverlaps(clusters_ranges, siW, maxgap = 100, ignore.strand = F)

##append the overlap info
clusters[queryHits(clus_olap)[!duplicated(queryHits(clus_olap))],]

toappend = as.data.frame(mcols(siW[subjectHits(clus_olap)[!duplicated(queryHits(clus_olap))]])) %>% 
  mutate(summ = paste(name, FDR, id, sep = ";"))

clusters$overlap = NA
clusters[queryHits(clus_olap)[!duplicated(queryHits(clus_olap))],]$overlap = toappend$summ

clusters = as_tibble(merge(clusters, mean_scores))

###count at the cluster level

regions_plus = promoters(clusters_ranges[strand(clusters_ranges)=="+"], upstream = 0,
                         downstream = 1000)
regions_minus = promoters(clusters_ranges[strand(clusters_ranges)=="-"], upstream = 0,
                          downstream = 1000)

#count for and rev separately - count ONLY RNA SEQ BAMS, check param
regcounts_plus <- regionCounts(bams_f, regions_plus, param=param)
regcounts_minus <- regionCounts(bams_r, regions_minus, param=param)

plusmat = sweep(x = assay(regcounts_plus), MARGIN = 2, STATS = regcounts_plus$totals, FUN = "/")
minmat = sweep(x = assay(regcounts_minus), MARGIN = 2, STATS = regcounts_minus$totals, FUN = "/")

regions_plus$ctlcounts = rowMeans(plusmat[,c(1:3)])
regions_plus$auxcounts = rowMeans(plusmat[,c(4:6)])

regions_minus$ctlcounts = rowMeans(minmat[,c(1:3)])
regions_minus$auxcounts = rowMeans(minmat[,c(4:6)])

regions_all = c(regions_plus, regions_minus)

clusters = as_tibble(merge(clusters, as.data.frame(regions_all)[,c(7,12,13)]))

clusters_final = clusters[order(clusters$cluster_mean_score, decreasing = T),] %>% 
  mutate(t1 = log2(ctlcounts*1e6), t2 = log2(auxcounts*1e6))

saveRDS(clusters_final, file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/ZC3H4_final/clusters_final.RDS")

###### affected by zc3h4
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

clusters_final$quartile = cut(clusters_final$cluster_mean_score, 
                              breaks = rev(quantile(clusters_final$cluster_mean_score, probs = seq(0,1,0.20))), 
                              labels = c(1,2,3,4,5), 
                              include.lowest = T)

clusters_final$quartile = factor(clusters_final$quartile, levels = c(5,4,3,2,1))

clusters_final[,c(13,14,17)]  %>% 
  pivot_longer(., cols = ends_with("counts"), names_repair = "unique") %>%
  mutate_at("name", as.factor) %>%
  mutate(name = fct_relevel(name, "ctlcounts", "auxcounts")) %>%
  ggplot(aes(quartile, log2(value*1e6))) + 
  geom_violin(aes(fill = name), width=0.6, alpha = 0.4) + 
  geom_boxplot(aes(fill = name), width=0.1, color="grey", alpha=0.8, position=position_dodge(0.6)) + 
  theme_classic() + 
  scale_fill_manual(values=cbbPalette) + 
  xlab("")

unlist(lapply(5:1, function(x){
  wilcox.test(log2(clusters_final[clusters_final$quartile == x,]$ctlcounts+1), 
              log2(clusters_final[clusters_final$quartile == x,]$auxcounts+1))$p.val
}))

#### consolidated peaks, individual

peaks_table = as.data.frame(peaks)

peaks_olap = findOverlaps(peaks, siW, maxgap = 100, ignore.strand = F)

##append the overlap info
peaks_table[queryHits(peaks_olap)[!duplicated(queryHits(peaks_olap))],]

#toappend = as.data.frame(mcols(siW[subjectHits(peaks_olap)[!duplicated(queryHits(peaks_olap))]])) %>% 
#  mutate(summ = paste(name, FDR, id, sep = ";"))

toappend = as.data.frame(mcols(siW[subjectHits(peaks_olap)[!duplicated(queryHits(peaks_olap))]])) %>% 
  mutate(summ = paste(name, FDR, id, sep = ";"))

peaks_table$overlap = NA
peaks_table[queryHits(peaks_olap)[!duplicated(queryHits(peaks_olap))],]$overlap = toappend$summ

peaks_table = as_tibble(merge(peaks_table, mean_scores))

peaks_table = peaks_table[order(peaks_table$adj_score, decreasing = T),]

saveRDS(peaks_table, 
        file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/ZC3H4_final/idr_peaks_wolap.RDS")

#make siW w eCLIP
peaks_olap = findOverlaps(siW, peaks, maxgap = 100, ignore.strand = F)
mcols_to_add = mcols(peaks[(as.data.frame(peaks_olap) %>% group_by(queryHits))$subjectHits])

#to_add = cbind(as.data.frame(peaks_olap), mcols_to_add) %>%
#  group_by(queryHits) %>%
#  mutate(eCLIP_peak_id = paste0("id", subjectHits-1)) %>%
#  dplyr::select(-c(log10pval, log2FC, rpkm, subjectHits)) %>%
#  summarise_all(paste, collapse = ",") %>%
#  mutate(siW_id = paste0("id", queryHits))

to_add = cbind(as.data.frame(peaks_olap), mcols_to_add) %>%
  group_by(queryHits) %>%
  mutate(eCLIP_peak_id = paste0("id", subjectHits)) %>%
  dplyr::select(-c(log10pval, log2FC, rpkm, subjectHits)) %>%
  summarise_all(paste, collapse = ",") %>%
  mutate(siW_id = paste0("id", queryHits))

final_siW = makeGRangesFromDataFrame(merge(as.data.frame(siW), to_add, by.x = 'id', by.y = 'siW_id', all.x = T) %>%
  dplyr::select(-c(queryHits)), keep.extra.columns = T)

final_siW = final_siW[order(final_siW$adj_FDR, decreasing = T)]

## counts
#siw_counts = summarizeOverlaps(features = final_siW,
#                               reads = c("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/07/07.bam",
#                                         "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/15/15.bam",
#                                         "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/23/23.bam"), 
#                               inter.feature = FALSE,
#                               singleEnd=FALSE, 
#                               fragments=F, ignore.strand = F, 
#                               BPPARAM = MulticoreParam(),
#                               strandMode = 2)
saveRDS(siw_counts, "ZC3H4_final/siw_counts_WT.RDS")
rpkms = rowMeans(sweep(sweep(assay(siw_counts), MARGIN = 2, FUN = "/", STATS = c(55527162, 51626298, 58757683)), 
              MARGIN = 1, FUN = "/", STATS = width(rowRanges(siw_counts)))*1e9)

final_siW$WT_rpkm = rpkms

saveRDS(final_siW, file = "ZC3H4_final/siW_toaddtoCLIPxl.RDS")


#####rep1 and rep2, in case you are asked

#stem = "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/clipper/using_siW_w_coding/nomismatch/outputs/21Jun/"
stem = "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/clipper/10Mqy_ttome/post_reseq/outputs/26Oct/"

rep1info = read.table(paste0(stem, "rep1_sig.bed.full"), sep = "\t", header = F)
rep2info = read.table(paste0(stem, "rep2_sig.bed.full"), sep = "\t", header = F)

rep1ranges = customimportbed6(paste0(stem, "rep1_sig.bed"))
rep2ranges = customimportbed6(paste0(stem, "rep2_sig.bed"))

###overlap with siW - rep1

rep1_olap = findOverlaps(rep1ranges, siW, maxgap = 100, ignore.strand = F)

##append the overlap info
rep1info[queryHits(rep1_olap)[!duplicated(queryHits(rep1_olap))],]

# if using Viviana's list
#toappend = as.data.frame(mcols(siW[subjectHits(rep1_olap)[!duplicated(queryHits(rep1_olap))]])) %>% 
#  mutate(summ = paste(V4, V7, V8, V9, sep = ";"))

toappend = as.data.frame(mcols(siW[subjectHits(rep1_olap)[!duplicated(queryHits(rep1_olap))]])) %>% 
  mutate(summ = paste(name, FDR, id, sep = ";"))

rep1info$overlap = NA
rep1info[queryHits(rep1_olap)[!duplicated(queryHits(rep1_olap))],]$overlap = toappend$summ

rep1info_final = as.data.frame(rep1info) %>% 
  dplyr::select(-c("V4", "V7", "V8", "V9", "V10"))

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
saveRDS(rep1info_final, file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_26Oct/rep1info_final_28Aug.RDS")

#### rep 2

rep2_olap = findOverlaps(rep2ranges, siW, maxgap = 100, ignore.strand = F)

##append the overlap info
rep2info[queryHits(rep2_olap)[!duplicated(queryHits(rep2_olap))],]

# if using Viviana's list
#toappend = as.data.frame(mcols(siW[subjectHits(rep2_olap)[!duplicated(queryHits(rep2_olap))]])) %>% 
#  mutate(summ = paste(V4, V7, V8, V9, sep = ";"))

toappend = as.data.frame(mcols(siW[subjectHits(rep2_olap)[!duplicated(queryHits(rep2_olap))]])) %>% 
  mutate(summ = paste(name, FDR, id, sep = ";"))

rep2info$overlap = NA
rep2info[queryHits(rep2_olap)[!duplicated(queryHits(rep2_olap))],]$overlap = toappend$summ

rep2info_final = as.data.frame(rep2info) %>% 
  dplyr::select(-c("V4", "V7", "V8", "V9", "V10"))

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
saveRDS(rep2info_final, file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_26Oct/rep2info_final_28Aug.RDS")


######in case you need to come back...
save.image("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/ZC3H4_final/peaks_info_final.RData")


########### FIg 2F

as.data.frame(final_siW) %>%
  dplyr::select(WT_rpkm, eCLIP_peak_id) %>%
  mutate(RPKM2 = case_when(is.na(WT_rpkm) ~ 0.001,
                           .default = WT_rpkm)) %>%
  mutate(olap = !is.na(eCLIP_peak_id)) %>%
  ggplot(aes(x = olap, y = log2(RPKM2))) + 
  geom_violin() + 
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  theme_classic()

metadata_df = as.data.frame(final_siW) %>%
  dplyr::select(WT_rpkm, eCLIP_peak_id) %>%
  mutate(RPKM2 = case_when(is.na(WT_rpkm) ~ 0,
                           .default = WT_rpkm)) %>%
  mutate(olap = !is.na(eCLIP_peak_id)) 

wilcox.test(metadata_df[metadata_df$olap == T,]$RPKM2,
            metadata_df[metadata_df$olap == F,]$RPKM2)$p.val


as.data.frame(final_siW) %>%
  dplyr::select(WT_rpkm, eCLIP_peak_id) %>%
  mutate(RPKM2 = case_when(is.na(WT_rpkm) ~ 0.001,
                           .default = WT_rpkm)) %>%
  mutate(olap = !is.na(eCLIP_peak_id)) %>%
  ggplot(aes(x = olap, y = log2(RPKM2+1))) + 
  geom_boxplot(width = 0.5, outlier.shape = NA, notch = T,
               notchwidth = 0.6, staplewidth = 0.5) + 
  theme_classic() + ylim(c(0,1))



###### 22May 2025
### clusters with no overlap with siW - plot those

clusters_noolap = clusters_final[is.na(clusters_final$overlap),]

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

clusters_noolap$quartile = cut(clusters_noolap$cluster_mean_score, 
                              breaks = rev(quantile(clusters_noolap$cluster_mean_score, probs = seq(0,1,0.20))), 
                              labels = c(1,2,3,4,5), 
                              include.lowest = T)

clusters_noolap$quartile = factor(clusters_noolap$quartile, levels = c(5,4,3,2,1))

clusters_noolap[,c(13,14,17)]  %>% 
  pivot_longer(., cols = ends_with("counts"), names_repair = "unique") %>%
  mutate_at("name", as.factor) %>%
  mutate(name = fct_relevel(name, "ctlcounts", "auxcounts")) %>%
  ggplot(aes(quartile, log2(value*1e6))) + 
  geom_violin(aes(fill = name), width=0.6, alpha = 0.4) + 
  geom_boxplot(aes(fill = name), width=0.1, color="grey", alpha=0.8, position=position_dodge(0.6)) + 
  theme_classic() + 
  scale_fill_manual(values=cbbPalette) + 
  xlab("")

unlist(lapply(5:1, function(x){
  wilcox.test(log2(clusters_noolap[clusters_noolap$quartile == x,]$ctlcounts+1), 
              log2(clusters_noolap[clusters_noolap$quartile == x,]$auxcounts+1))$p.val
}))



















