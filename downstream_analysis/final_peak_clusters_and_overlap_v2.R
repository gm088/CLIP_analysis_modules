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
#siW = readRDS("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/T3_ghostbusters/ZC3H4_supp_tx_PROCAP_FANTOM.RDS")
siW = readRDS("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/ZC3H4_supp_tx_barebones.RDS")

names(mcols(siW)) = c("name", "FDR", "class", "WT_rpkm", "adj_FDR", "annot", "class", "lnc")
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


#######################
rel_siW = siW[subjectHits(clus_olap)[!duplicated(queryHits(clus_olap))]]
toappend = paste0(as.vector(seqnames(rel_siW)), ":", start(rel_siW), "-", end(rel_siW), ":", as.vector(strand(rel_siW)))

#toappend = as.data.frame(mcols(siW[subjectHits(clus_olap)[!duplicated(queryHits(clus_olap))]])) %>% 
#  mutate(summ = paste(name, FDR, id, sep = ";"))
#######################

clusters$overlap = NA
clusters[queryHits(clus_olap)[!duplicated(queryHits(clus_olap))],]$overlap = toappend

clusters = as_tibble(merge(clusters, mean_scores))

write.table(clusters[order(clusters$cluster_mean_score, decreasing = T),], 
            sep = "\t", file = "ZC3H4_final/v2/clusters_append.txt", 
            row.names = F, col.names = T, quote = F)


#### consolidated peaks, individual

peaks_table = as.data.frame(peaks)

peaks_olap = findOverlaps(peaks, siW, maxgap = 100, ignore.strand = F)

##append the overlap info
peaks_table[queryHits(peaks_olap)[!duplicated(queryHits(peaks_olap))],]

rel_siW2 = siW[subjectHits(peaks_olap)[!duplicated(queryHits(peaks_olap))]]
toappend = paste0(as.vector(seqnames(rel_siW2)), ":", start(rel_siW2), "-", end(rel_siW2), ":", as.vector(strand(rel_siW2)))

#######################

peaks_table$overlap = NA
peaks_table[queryHits(peaks_olap)[!duplicated(queryHits(peaks_olap))],]$overlap = toappend

peaks_table = as_tibble(merge(peaks_table, mean_scores))

peaks_table = peaks_table[order(peaks_table$adj_score, decreasing = T),]

write.table(peaks_table$overlap, "ZC3H4_final/v2/idr_peaks_append.txt", quote = F, row.names = F)

saveRDS(peaks_table, 
        file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/ZC3H4_final/v2/idr_peaks_wolap.RDS")

#make siW w eCLIP
peaks_olap = findOverlaps(siW, peaks, maxgap = 100, ignore.strand = F)

# I want to add the eCLIP peak id and the coordinates
peaks_to_append = peaks[(as.data.frame(peaks_olap) %>% group_by(queryHits))$subjectHits]
coords_peaks_to_append = paste0(as.vector(seqnames(peaks_to_append)), ":", start(peaks_to_append), "-", 
                                end(peaks_to_append), ":", as.vector(strand(peaks_to_append)))

mcols_to_add = cbind(names(peaks_to_append), coords_peaks_to_append)

to_add = cbind(as.data.frame(peaks_olap), mcols_to_add) %>%
  group_by(queryHits) %>%
  mutate(eCLIP_peak_id = V1) %>%
  dplyr::select(-c(V1, subjectHits)) %>%
  summarise_all(paste, collapse = ";") %>%
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
#saveRDS(siw_counts, "ZC3H4_final/siw_counts_WT.RDS")
#rpkms = rowMeans(sweep(sweep(assay(siw_counts), MARGIN = 2, FUN = "/", STATS = c(55527162, 51626298, 58757683)), 
#              MARGIN = 1, FUN = "/", STATS = width(rowRanges(siw_counts)))*1e9)
#
#final_siW$WT_rpkm = rpkms

saveRDS(final_siW, file = "ZC3H4_final/v2/siW_table_s2.RDS")
write.table(as.data.frame(final_siW), file = "ZC3H4_final/v2/siw.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)



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



















