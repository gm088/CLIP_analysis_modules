source("feat_extr_kmer_functions.R")
source("functions_2.R")
source("functions.R")
library(GenomicAlignments)
library(BiocParallel)

###### metaplot
BList = c(BList, customimportbed6("/hpcnfs/data/GN2/gmandana/annotation/small_RNA_tRNA_RPL_RPPH_bed6.bed"))

bwstem = "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/"

# import siW
bw_f = paste0(bwstem, c("13_for.bw", "94_merged_for.bw", "95_merged_for.bw"))
bw_r = paste0(bwstem, c("13_rev.bw", "94_merged_rev.bw", "95_merged_rev.bw"))


siW = readRDS("ZC3H4_final/v2/siW_table_s2.RDS")
reltx = siW[!is.na(siW$eCLIP_peak_id)]
notreltx = subsetByOverlaps(siW[is.na(siW$eCLIP_peak_id)], BList, invert = T)

prof_siW = covmat_noscalereg_w_flanks(bwfiles = list(bw_f, bw_r),
                                      regions = reltx[width(reltx) > 1000],
                                      numbins = 1000, 
                                      ups = 1000, downs = 1000,
                                      upbins = 100, downbins = 100)

prof_siW_noCLIP = covmat_noscalereg_w_flanks(bwfiles = list(bw_f, bw_r),
                                      regions = notreltx[width(notreltx) > 1000],
                                      numbins = 1000, 
                                      ups = 1000, downs = 1000,
                                      upbins = 100, downbins = 100)

### metaplotter
colors = c("#522302", "#730a0a", "#adadad","#000000")

#labels = gsub(rownames(prof_siW), pattern = "_for", replacement = "")
data_prof = metaplotter(matrix = prof_siW, metadata = c("CLIP_pos_INPUT", "CLIP_pos_ZC3H4_R1", "CLIP_pos_ZC3H4_R2"), return_df = T)
data_prof2 = metaplotter(matrix = prof_siW_noCLIP, metadata = c("INPUT", "ZC3H4_R1", "ZC3H4_R2"), return_df = T)

rbind(data_prof, data_prof2) %>%
  ggplot(aes(x = bin, y = value, fill = group)) + 
  geom_line(aes(color = group)) + 
  geom_ribbon(aes(y = value, ymin = value+(value_sd), 
                  ymax = value-(value_sd), group = group), 
              alpha = 0.5) + 
  theme_classic() + xlab("") + ylab("") + 
  scale_x_continuous(breaks = c(100,1100), labels = c("TSS", "TES")) + 
  scale_color_manual(values = colors) + scale_fill_manual(values = colors)


##### coutning

bams = c("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/mapped/14/input2_clean_sorted.bam",
         "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/mapped/13/input1_clean_sorted.bam",
         "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/mapped/IP7_merged/94_merged.bam",
         "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/mapped/95/DCIP_clean_sorted.bam")

#siw_counts = summarizeOverlaps(features = reltx,
#                               reads = bams, 
#                               inter.feature = FALSE,
#                               singleEnd=FALSE, 
#                               fragments=F, ignore.strand = F, 
#                               BPPARAM = MulticoreParam(),
#                               strandMode = 2)
#
#siw_counts_noCLIP = summarizeOverlaps(features = notreltx,
#                                      reads = bams, 
#                                      inter.feature = FALSE,
#                                      singleEnd=FALSE, 
#                                      fragments=F, ignore.strand = F, 
#                                      BPPARAM = MulticoreParam(),
#                                      strandMode = 2)

#save.image("ZC3H4_final/no_overlap_metaplot/this.RData")

libsizes = c(19520190,17736868,22567492,21522129)

rpkms = sweep(sweep(assay(siw_counts), MARGIN = 2, FUN = "/", STATS = libsizes), 
              MARGIN = 1, FUN = "/", STATS = width(rowRanges(siw_counts)))*1e9

rpkms_noCLIP = sweep(sweep(assay(siw_counts_noCLIP), MARGIN = 2, FUN = "/", STATS = libsizes), 
                              MARGIN = 1, FUN = "/", STATS = width(rowRanges(siw_counts_noCLIP)))*1e9

avg_rpkm = cbind(rowMeans(rpkms[,c(1,2)]), rowMeans(rpkms[,c(3,4)]))
avg_rpkm_noCLIP = cbind(rowMeans(rpkms_noCLIP[,c(1,2)]), rowMeans(rpkms_noCLIP[,c(3,4)]))

boxplot(avg_rpkm_noCLIP, notch = T, outline = F)

wilcox.test(avg_rpkm[,1], avg_rpkm[,2], paired = T)$p.value
wilcox.test(avg_rpkm_noCLIP[,1], avg_rpkm_noCLIP[,2], paired = T)$p.value




### 397 discarded

avg_rpkm_noCLIP_nozeros = avg_rpkm_noCLIP[rowSums(avg_rpkm_noCLIP == 0) == 0,]

wilcox.test(avg_rpkm_noCLIP_nozeros[,1], avg_rpkm_noCLIP_nozeros[,2], paired = T)$p.value

wilcox.test(avg_rpkm[,1], avg_rpkm[,2], paired = T)$p.value

samples = sapply(1:1000, function(x){wilcox.test(sample(avg_rpkm_noCLIP_nozeros[,1], 1000), 
                                       sample(avg_rpkm_noCLIP_nozeros[,2], 1000), paired = T)$p.value})

samples2 = sapply(1:1000, function(x){wilcox.test(sample(avg_rpkm[,1], 200), 
                                                 sample(avg_rpkm[,2], 200), paired = T)$p.value})







