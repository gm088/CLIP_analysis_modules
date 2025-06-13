source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/Shannahan/functions.R")
library(RColorBrewer)

b_f_se = c("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/13_SE_fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/94_merged_read2_SE_fwd.bam",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/DRB_DMSO_SE_fwd.bam",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/DRB_5min_merged_SE_fwd.bam",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/CSTF_2_NT_R1_merged_SE_fwd.bam",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/CSTF2_AUX_merged_SE_fwd.bam",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/PNUTS_R1_merged_SE_fwd.bam",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/PNUTS_R2_merged_SE_fwd.bam")

b_r_se = c("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/13_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/94_merged_read2_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/DRB_DMSO_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/DRB_5min_merged_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/CSTF_2_NT_R1_merged_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/CSTF2_AUX_merged_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/PNUTS_R1_merged_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/PNUTS_R2_merged_SE_rev.bam")

win_width = 100
param = readParam(minq=1, pe="none")

#regions = import('/hpcnfs/data/GN2/gmandana/fc2_m1_rbm26peaks.bed')
#regions = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_26Oct/reproducible_peaks_reseq_rpkm.bed")
#names(regions) = paste0("id", 1:length(regions))

regions = regs_annot

regions_plus = regions[strand(regions)=="+"]
regions_minus = regions[strand(regions)=="-"]

#count for and rev separately
regcounts_plus <- regionCounts(b_f_se, regions_plus, param=param)
regcounts_minus <- regionCounts(b_r_se, regions_minus, param=param)

norm_plus = sweep(assay(regcounts_plus), 2, regcounts_plus$totals, FUN = "/") * 1e6
norm_minus = sweep(assay(regcounts_minus), 2, regcounts_minus$totals, FUN = "/") * 1e6

rpms = rbind(norm_plus, norm_minus)

colnames(rpms) = c("Input", "IP1", "DRB_DMSO", "DRB_5min", 
                     "CSTF2_NT", "CSTF2_AUX", "PNUTS_R1", "PNUTS_R2")

rownames(rpms) = c(names(regions_plus), names(regions_minus))

pheatmap(log2(rpms+1),
         cluster_cols = F,
         show_rownames = F)


#####one sec
regs = c(regions_plus, regions_minus)


pheatmap(log2(rpms[names(regs[regs$siW == T,]),c(1,2,4,8,5,6)]+1),
         cluster_cols = F,
         show_rownames = F)



###cluster and downstream count

####pivot
clusters = as.data.frame(regions) %>% 
  group_by(V4) %>% 
  mutate(cluster_start = min(start)) %>% 
  mutate(cluster_end = max(end)) %>% 
  mutate(coords = paste(seqnames, paste(start, end, sep = "-"), sep = ":")) %>% 
  dplyr::select(-c("start", "end", "width")) %>% 
  summarise_all(paste, collapse = ",")

clusters$seqnames = unlist(lapply(clusters$seqnames, 
                                  function(x){str_split(x, ",")[[1]][[1]]}))
clusters$strand = unlist(lapply(clusters$strand, 
                                function(x){str_split(x, ",")[[1]][[1]]}))
clusters$cluster_start = as.numeric(unlist(lapply(clusters$cluster_start, 
                                                  function(x){str_split(x, ",")[[1]][[1]]})))
clusters$cluster_end = as.numeric(unlist(lapply(clusters$cluster_end, 
                                                function(x){str_split(x, ",")[[1]][[1]]})))

clusters_ranges = makeGRangesFromDataFrame(clusters,
                                           start.field = "cluster_start",
                                           end.field = "cluster_end",
                                           strand.field = "strand",
                                           keep.extra.columns = T)

clusters_ranges_p = promoters(clusters_ranges[strand(clusters_ranges)=="+"],
                              upstream = 0, downstream = 1000)
clusters_ranges_m = promoters(clusters_ranges[strand(clusters_ranges)=="-"],
                              upstream = 0, downstream = 1000)

#count for and rev separately
regcounts_plus <- regionCounts(bams_f, clusters_ranges_p, param=param)
regcounts_minus <- regionCounts(bams_r, clusters_ranges_m, param=param)

norm_plus = sweep(assay(regcounts_plus), 2, regcounts_plus$totals, FUN = "/") * 1e6
norm_minus = sweep(assay(regcounts_minus), 2, regcounts_minus$totals, FUN = "/") * 1e6

plus_ctlcounts = rowMeans(norm_plus[,c(1:3)])
plus_auxcounts = rowMeans(norm_plus[,c(4:6)])

minus_ctlcounts = rowMeans(norm_minus[,c(1:3)])
minus_auxcounts = rowMeans(norm_minus[,c(4:6)])

totmat = rbind(cbind(plus_ctlcounts, plus_auxcounts), 
      cbind(minus_ctlcounts, minus_auxcounts))

boxplot(log2(totmat[regs$class == "intron",]+1))

rpms = rbind(norm_plus, norm_minus)

colnames(rpms) = c("Input", "IP1", "DRB_DMSO", "DRB_5min", 
                   "CSTF2_NT", "CSTF2_AUX", "PNUTS_R1", "PNUTS_R2")

rownames(rpms) = c(names(clusters_ranges_p), names(clusters_ranges_m))

#####one sec
regs = c(clusters_ranges_p, clusters_ranges_m)

pheatmap(log2(rpms[names(regs[regs$class == "3.UTR",]),c(1,2,4,8,5,6)]+1),
         cluster_cols = F,
         show_rownames = F)





