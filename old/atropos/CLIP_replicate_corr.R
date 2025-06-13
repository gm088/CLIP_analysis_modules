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
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = loadDb("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/txdb.sqlite"); seqlevels(txdb) = paste0("chr", seqlevels(txdb))

alltx = transcripts(txdb, use.names = T, columns=c("TXNAME", "TXTYPE"))
coding = alltx[alltx$TXTYPE == "protein_coding"]

bams_f = c("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/13_SE_fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/14_SE_fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/94_SE_fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/95_SE_fwd.bam")

bams_r = c("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/13_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/14_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/94_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/SE/95_SE_rev.bam")

samples = 1:4
Desc = c("INPUT","INPUT", "IP", "IP")
metadata = data.frame(samples, Desc)
desmat <- model.matrix(~factor(metadata$Desc, levels = c("INPUT", "IP")))
colnames(desmat) = c("INPUT", "IP")
param = readParam(minq=10, pe="none")

win_width = 500

binned_f <- windowCounts(bams_f, bin=TRUE, width=win_width, param=param, BPPARAM=MulticoreParam(), filter=100)
binned_r <- windowCounts(bams_r, bin=TRUE, width=win_width, param=param, BPPARAM=MulticoreParam(), filter=100)

allcounts = rbind(assay(binned_f), assay(binned_r))
colnames(allcounts) = c("INPUT-1","INPUT-2", "IP-1", "IP-2")

libsizes = regcounts_plus$totals + regcounts_minus$totals

normbinned = sweep(sweep(allcounts, 2, libsizes*1e-6, FUN = "/"), 1, reg_widths, FUN = "*")
#normbinned = normbinned[rowSums(normbinned) < quantile(rowSums(normbinned), probs = seq(0,1,0.01))[99],]

data.pca = princomp(cor(normbinned))
fviz_eig(data.pca, addlabels = TRUE)
fviz_pca_var(data.pca, col.var = "black")

###counting in predefined regions
## which regions to use? - genes, zc3h4 repressed tx, or peaks?

source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")
peaks_file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/files_for_svm/new_top_peaks_clipper.bed"
unaff_file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/files_for_svm/all5prime_expressed_unique_500bp_notdiffexpr.bed"
siW = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/new_siW/all_siW_new.bed"
unaff_pcg = subsetByOverlaps(import(unaff_file), coding)
siW_regions = customimportbed6(siW)
peak_regions = customimportbed6(peaks_file)

regions = peak_regions
regions_plus = regions[strand(regions)=="+"]
regions_minus = regions[strand(regions)=="-"]

regcounts_plus <- regionCounts(bams_f, regions_plus, param=param)
regcounts_minus <- regionCounts(bams_r, regions_minus, param=param)

allcounts = rbind(assay(regcounts_plus), assay(regcounts_minus))
colnames(allcounts) = c("INPUT_1","INPUT_2", "IP_1", "IP_2")

libsizes = regcounts_plus$totals + regcounts_minus$totals
reg_widths = c(width(regions_plus), width(regions_minus)) 

normbinned = sweep(sweep(allcounts, 2, libsizes*1e-9, FUN = "/"), 1, reg_widths, FUN = "/")
normbinned = normbinned[rowSums(normbinned) < quantile(rowSums(normbinned), probs = seq(0,1,0.01))[99],]

data.pca = princomp(cor(normbinned))
fviz_eig(data.pca, addlabels = TRUE)
fviz_pca_var(data.pca, col.var = "black")

corr_fn(normbinned_mat = normbinned, "peak_regions_correlations.pdf")

corr_fn = function(normbinned_mat, outname){
  
  
  pos = log10(median(rowSums(normbinned_mat)))
  pdf(outname, width = 10, height = 5)
  
  correlation_1 = cor(log10(normbinned_mat[,3] + 1), log10(normbinned_mat[,4] + 1))
  #correlation scatter
  plot1 = as.data.frame(normbinned_mat) %>% 
    ggplot(aes(log10(IP_1), log10(IP_2))) + 
    geom_point()  + theme_classic() + 
    xlab("log10(RPKM IP1)") + ylab("log10(RPKM IP2)") + 
    annotate("text", x=0, y=2, label= bquote(R^2 ~ "=" ~ .(correlation_1^2)))
  print(plot1)
  
  correlation_2 = cor(log10(normbinned_mat[,1] + 1), log10(normbinned_mat[,2] + 1))
  
  plot2 = as.data.frame(normbinned_mat) %>% 
    ggplot(aes(log10(INPUT_1), log10(INPUT_2))) + 
    geom_point()  + theme_classic() + 
    xlab("log10(RPKM INPUT1)") + ylab("log10(RPKM INPUT2)") + 
    annotate("text", x=0, y=2, label= bquote(R^2 ~ "=" ~ .(correlation_2^2)))
  print(plot2)
  
  data.pca = princomp(cor(normbinned_mat))
  #fviz_eig(data.pca, addlabels = TRUE)
  plot3 = fviz_pca_var(data.pca, col.var = "black")
  print(plot3)
  dev.off()
  
}




