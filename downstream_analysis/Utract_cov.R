library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
#library('biomaRt')
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/Shannahan/functions.R")

bedfile = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_26Oct/reproducible_peaks_reseq_rpkm.bed"
clipper_peaks = customimportbed6(bedfile)
mysiW = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/T3_ghostbusters/siW__better_annotated_20Sep.bed")

#I want the procap, IP7, PROSEQ, TT-seq
b_f = c("/hpcnfs/data/GN2/gmandana/take2project/PROSEQ/bigwigs/28_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/bigwigs_RPKM/23_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/bigwigs_RPKM/24_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/13_r2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/94_merged_read2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/95_r2_for.bw",
        #"/hpcnfs/data/GN2/gmandana/GSM4593586_L_EGFP_rep1_tt_corr_ff_noJncReads_plus.bw",
        #"/hpcnfs/data/GN2/gmandana/GSM4593592_L_RRP40_rep1_tt_corr_ff_noJncReads_plus.bw",
        #"/hpcnfs/data/GN2/gmandana/take2project/LisCap/bigwigs/SRR22522159_fwd_5p.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/DCIP_w_multimap_for.bw")

b_r = c("/hpcnfs/data/GN2/gmandana/take2project/PROSEQ/bigwigs/28_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/bigwigs_RPKM/23_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/bigwigs_RPKM/24_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/13_r2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/94_merged_read2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/95_r2_rev.bw",
        #"/hpcnfs/data/GN2/gmandana/GSM4593586_L_EGFP_rep1_tt_corr_ff_noJncReads_minus.bw",
        #"/hpcnfs/data/GN2/gmandana/GSM4593592_L_RRP40_rep1_tt_corr_ff_noJncReads_minus.bw",
        #"/hpcnfs/data/GN2/gmandana/take2project/LisCap/bigwigs/SRR22522159_rev_5p.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/DCIP_w_multimap_rev.bw")

#tags = c("PROSEQ", "IP7", "DCIP", "TTseq_ctl", "TTseq_siRRP40", "PROCAP", "DCIP_w_multimapp")
tags = c("PROSEQ", "4Su_ctl", "4Su_AUX", "INPUT", "IP7", "DCIP", "DCIP_w_multimapp")
cbbPalette <- c("#000000", "#D19696", "#D80000", "#C2C2C2", 
                "#0B9200", "#0B9200", "#007A92")

#pU = customimportbed6("U_tracts/clust_Utract_hits.bed")
pU = customimportbed6("U_tracts/bg3_Utract_hits.bed")
#pU = customimportbed6("U_tracts/randomgene_hits2.bed")

#these are your ranges to plot at cluster level

pU_aligned_prof = getcovmat2(list(b_f, b_r), 
                                promoters(pU, upstream = 500, downstream = 500), 
                                nbins = 1000, norm = F)

#pU_aligned_prof2 = getcovmat(list(b_f, b_r), 
#                            promoters(pU, upstream = 500, downstream = 500), 
#                            nbins = 1000, norm = T)

#rownames(pU_aligned_prof) = tags

melt(pU_aligned_prof) %>%
  ggplot(aes(x = Var2, y = value, group = Var1)) + 
  geom_line(aes(color = Var1)) + 
  scale_color_manual(values=cbbPalette) + 
  theme_classic() + xlab("") + ylab("")


## small sanity check...
profs_ip7_f = bigwig_collect(b_f[[2]], promoters(pU[strand(pU) == "+"], upstream = 2000, downstream = 2000),
                             norm = F, nbins = 400)

### remove anything within coding gene boundaries...
txdb = loadDb("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/txdb.sqlite")
seqlevels(txdb) = paste0("chr", seqlevels(txdb))
alltx = transcripts(txdb, use.names = T, columns = c('TXTYPE', 'GENEID'))
pctx = alltx[alltx$TXTYPE == 'protein_coding']

proseq_peaks = customimportbed6("/hpcnfs/data/GN2/gmandana/take2project/PROSEQ/macs/round3/MACS_peaks.narrowPeak")
























