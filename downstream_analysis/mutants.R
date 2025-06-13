source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/Shannahan/functions.R")
library(RColorBrewer)


BL = import("/hpcnfs/data/GN2/gmandana/annotation/small_RNA_tRNA_RPL_RPPH_bed6.bed")
####ranges to consider

mysiW = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/T3_ghostbusters/nctx_refinement/siW_refined_annotated_fantom.bed")
names(mcols(mysiW)) = c("name", "FDR", "class", "ctl_rpkm", "fantom", "class2")

eCLIP_peaks = subsetByOverlaps(customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_26Oct/reproducible_peaks_reseq_rpkm.bed"), 
                               BL, invert = T)

mysiW = mysiW[!(names(mysiW) %in% names(subsetByOverlaps(mysiW, BL)[width(subsetByOverlaps(mysiW, BL)) < 1000]))]

opregs = readRDS("~/ENHANCEDCLIP/ZC3H4_final/opregs.RDS")

gencode = customimportbed6("/hpcnfs/data/GN2/gmandana/annotation/gencode.v39.bed")
SNHG_hist = gencode[grepl("SNHG|^H[0-9].*", gencode$V13)]

queryregs = subsetByOverlaps(opregs, SNHG_hist, invert = T, maxgap = 500)[1:1001]
#######

# metaplots

#######
bwdir="/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/"

b_f = paste0(bwdir, c("94_merged_for.bw", "delRGG_INPUT_R2_R2_for.bw", "13_for.bw",
        "WT_22May_R2_R2_for.bw", "WT_22May_R4_R2_for.bw",
        "delCCCH3_R2_R2_for.bw","delCCCH3_R3_R2_for.bw",
        "delRGG_R2_R2_for.bw", "delRGG_R1_R2_for.bw"))

b_r = paste0(bwdir, c("94_merged_rev.bw", "delRGG_INPUT_R2_R2_rev.bw", "13_rev.bw",
        "WT_22May_R2_R2_rev.bw", "WT_22May_R4_R2_rev.bw",
        "delCCCH3_R2_R2_rev.bw","delCCCH3_R3_R2_rev.bw",
        "delRGG_R2_R2_rev.bw", "delRGG_R1_R2_rev.bw"))

#whichbw = c(4,3,5,6,7,8,9,10,11)
whichbw = 1:length(b_r)

b_f_sel = b_f[whichbw]
b_r_sel = b_r[whichbw]

regions = queryregs

prof2 = getcovmat(list(b_f_sel, b_r_sel), 
                  promoters(regions, upstream = 100, downstream = 1000), 
                  nbins = 1100, norm = F)

saveRDS(prof2, "mutants/mutants_siW_prof_10May25.RDS")

labels = c("ZC3H4_R1", "mutants_INPUT", "ZC3H4_INPUT", "ZC3H4_WT_R1", "ZC3H4_WT_R2",
           "ZC3H4_delCCCH3_R1", "ZC3H4_delCCCH3_R2", "ZC3H4_delRGG_R1", "ZC3H4_delRGG_R2")

cbbPalette <- c("#85929E", "#17202A", "#82E0AA", "#186A3B", 
                "#EDBB99", "#F39C12", "#eb4034", "#5B2C6F", "#920000")

data_prof = metaplotter(prof2[2:9,], metadata = labels[2:9])

data_prof + scale_x_continuous(breaks = c(100,2100), labels = c("Peak Region Start", "2kb")) + 
  scale_color_manual(values = cbbPalette) + scale_fill_manual(values = cbbPalette)



######### coutnts for excel
bamdir="/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/mapped/"
bamdir2="/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/mapped/mutants_merged/"

bams = c(paste0(bamdir, "13/input1_clean_sorted.bam"),
         paste0(bamdir2, "WT_R2_merged.bam"),
         paste0(bamdir2, "WT_R4_merged.bam"),
         paste0(bamdir2, "delCCCH3_R2_merged.bam"),
         paste0(bamdir, "delCCCH3_R3/delCCCH3_R3.bam"),
         paste0(bamdir2, "delRGG_R2_merged.bam"),
         paste0(bamdir, "delRGG_R1/delRGG_R1.bam"))

#regs_counts = summarizeOverlaps(features = queryregs[1001] , 
#                                reads = bams, 
#                                inter.feature = FALSE,
#                                singleEnd=FALSE, 
#                                fragments=F, ignore.strand = F, 
#                                BPPARAM = MulticoreParam(),
#                                strandMode = 2)
#saveRDS(regs_counts, "mutants/top1000_regs_counts.RDS")

libsizes = c(35473736, 31356908, 27715846, 
             22751526, 18249712, 23391426, 32500630)

rpkms = sweep(sweep(assay(regs_counts), 
                    MARGIN = 1, STATS = width(rowRanges(regs_counts)), FUN = "/"),
              MARGIN = 2, STATS = libsizes, FUN = "/")*1e9

final_rpkm = cbind(rpkms[,1], 
                   rowMeans(rpkms[,c(2,3)]), 
                   rowMeans(rpkms[,c(4,5)]),
                   rowMeans(rpkms[,c(6,7)]))

boxplot(log2(final_rpkm+1), outline = F, notch = T)
wilcox.test(log2(final_rpkm[,4]+1), log2(final_rpkm[,2]+1))

colnames(final_rpkm) = c("INPUT", "WT", "delCCCH3", "delRGG")

total_df = cbind(as.data.frame(regions)[,c(1,2,3,4,5)], final_rpkm)

write.table(total_df, file = "mutants/table_s5_1.txt", quote = F, sep = "\t",
            row.names = F, col.names = T)

















