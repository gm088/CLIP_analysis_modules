source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/Shannahan/functions.R")
library(RColorBrewer)

b_f = c("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/bigwigs_RPKM/noncod/23_noncod_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/bigwigs_RPKM/noncod/24_noncod_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/PROSEQ/bigwigs/28_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/95_r2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/13_r2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/DRB_DMSO_r2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/DRB_0min_r2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/DRB_5min_merged_read2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/CSTF_2_NT_R1_merged_read2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/CSTF2_AUX_merged_read2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/PNUTS_R1_merged_read2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/PNUTS_R2_merged_read2_for.bw")

b_r = c("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/bigwigs_RPKM/noncod/23_noncod_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/bigwigs_RPKM/noncod/24_noncod_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/PROSEQ/bigwigs/28_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/95_r2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/13_r2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/DRB_DMSO_r2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/DRB_0min_r2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/DRB_5min_merged_read2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/CSTF_2_NT_R1_merged_read2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/CSTF2_AUX_merged_read2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/PNUTS_R1_merged_read2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/PNUTS_R2_merged_read2_rev.bw")

BL = import("/hpcnfs/data/GN2/gmandana/annotation/small_RNA_tothrow_and_just_tRNA_hg38_encode_BL_bed4.bed")
####ranges to consider

mysiW = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/T3_ghostbusters/siW__better_annotated_20Sep.bed")
mysiW = subsetByOverlaps(mysiW, BL, invert = T)
names(mcols(mysiW)) = c("name", "ctl_rpkm", "pval", "FDR", "logFC")
names(mysiW) = paste0("id", 1:length(mysiW))
mysiW = mysiW[width(mysiW) > 500]

r_p = promoters(mysiW[strand(mysiW) == "+"], upstream = 1000, downstream = 20000)
r_m = promoters(mysiW[strand(mysiW) == "-"], upstream = 1000, downstream = 20000)

indices = c(5,4,6,8:12)

plusmat = do.call('rbind', lapply(indices, function(x){
  colMeans(bigwig_collect(b_f[x], regions = r_p, norm = F, nbins = 400))
}))

minmat = do.call('rbind', lapply(indices, function(x){
  colMeans(bigwig_collect(b_r[x], regions = r_m, norm = F, nbins = 400))
}))

#because they are already averaged over regions
covmat = (minmat + plusmat)/2

smoothed = do.call('rbind', lapply(1:nrow(covmat), function(x){
  unlist(slide(.x = covmat[x,], .f = mean, .before = 5, .after = 5))
}))

rownames(covmat) = c("Input", "IP2", "DRB_DMSO", "DRB_5min", 
                     "CSTF2_NT", "CSTF2_AUX", "PNUTS_R1", "PNUTS_R2")

matplot(covmat[2,], type = "l", ylim = c(0,800))
matplot(covmat[1,], type = "l", add = T)
matplot(covmat[3,], type = "l", add = T)
matplot(covmat[4,], type = "l", add = T)
matplot(covmat[5,], type = "l", add = T)
matplot(covmat[6,], type = "l", add = T)

matplot(t(covmat), type = "l")
legend(300, 900, legend = rownames(covmat), col = 1:6, lty = 1:6,
       cex=0.5, pch=1, pt.cex = 1)

matplot(covmat[3,], type = "l")






