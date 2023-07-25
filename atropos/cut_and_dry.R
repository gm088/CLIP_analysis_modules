source("/Users/IEO5559/enhancedClip/Shannahan/feat_extr_kmer_functions.R")
source("~/enhancedClip/functions_2.R")
source("/Users/IEO5559/enhancedClip/Shannahan/functions.R")
options(meme_bin = "/Users/IEO5559/bin/bin")

#import the files you need for overlapping
my_siW = customimportbed6("/Users/IEO5559/prog/bedfiles/mysiW/fixed_filter50/mysiW_fdr0.0001_gt1kb.bed")
clipper_peaks = customimportbed6("/Users/IEO5559/enhancedClip/peak_files_24Jul/clipper_new_rpkm_noucscrepeats_nosmallRNA.bed")
colnames(mcols(clipper_peaks)) = c("txname",
                           "adj_score",
                           "log10pval",
                           "log2FC",
                           "rpkm")
new_siW = customimportbed6("/Users/IEO5559/rbns/all_siW_new.bed")

#get the clipper peaks
CLIP_TUs = CLIP_systemTUs_get(clip_peaks = clipper_peaks, systemtus = systemTUs)

#background region selection
bg1 = systemTUs[systemTUs$FDR > 0.1]
bg2 = systemTUs[systemTUs$FDR > 0.1 & systemTUs$iscoding == T]
#divergent units not sensitive
bg3 = systemTUs[div[div$FDR > 0.1 & abs(div$logFC) < 0.5]$NDR_tag]
#enhancer units not sensitive
bg4 = systemTUs[enh[enh$FDR > 0.1 & abs(enh$logFC) < 0.5]$NDR_tag]
#divergent mates not sensitive
div_sensitive = subsetByOverlaps(div, new_siW)
div_sensitive_mates = div[div$NDR %in% div_sensitive$NDR & !(div$NDR_tag %in% div_sensitive$NDR_tag)]

bg5 = systemTUs[div_sensitive_mates$NDR_tag]; bg5 = bg5[bg5$FDR > 0.1]

print("###")

ms_pro = streme_wrapper(promoters(CLIP_TUs, upstream = 1, downstream = 200), 
                     promoters(bg5, upstream = 1, downstream = 200), 
                     write_dir = "Shannahan/test2",
                     order = 4,
                     minwidth = 4,
                     maxwidth = 10)


#don't alter this for now
ms_test = streme_wrapper(promoters(get_uniq_TSS(systemTUs[div_sensitive$NDR_tag], mode = "main"), upstream = 1, downstream = 200), 
                         promoters(get_uniq_TSS(bg5, mode = "main"), upstream = 1, downstream = 200), 
                         write_dir = "Shannahan/test",
                         order = 4,
                         minwidth = 4,
                         maxwidth = 30)

### run the searches
TSS_and_motif(outdir = "outputs/bg2", filename = "CLIP_TU_vs_bg2.pdf", 
              fg = CLIP_TUs, bg = bg2, tss_mode = "main", 
              width = 100, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4)

TSS_and_motif(outdir = "outputs/bg3", filename = "CLIP_TU_vs_bg3.pdf", 
              fg = CLIP_TUs, bg = bg3, tss_mode = "main", 
              width = 100, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4)

TSS_and_motif(outdir = "outputs/tophalf_vs_bg3", filename = "tophalf_vs_bg3.pdf", 
              fg = CLIP_TUs[1:400], bg = bg3, tss_mode = "main", 
              width = 100, shift_val = -2, minw = 4, maxw = 20, streme_order = 1:4)

TSS_and_motif(outdir = "outputs/bg5", filename = "CLIP_TU_vs_bg5.pdf", 
              fg = CLIP_TUs, bg = bg5, tss_mode = "main", 
              width = 100, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4)

TSS_and_motif(outdir = "outputs/bg5_pro", filename = "CLIP_TU_vs_bg5.pdf", 
              fg = CLIP_TUs, bg = bg5, tss_mode = "main", usePROSEQ = T,
              width = 100, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4)

TSS_and_motif(outdir = "outputs/tophalf_vs_bg3_pro", filename = "tophalf_vs_bg3_pro.pdf", 
              fg = CLIP_TUs[1:400], bg = bg3, tss_mode = "main", usePROSEQ = T,
              width = 100, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4)

TSS_and_motif(outdir = "outputs/div_sensitive", filename = "div_vs_div_sensitive.pdf", 
              fg = systemTUs[div_sensitive$NDR_tag], bg = bg5, tss_mode = "main", 
              width = 100, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4)

TSS_and_motif(outdir = "outputs/tophalf_vs_bg5", filename = "tophalf_vs_bg5.pdf", 
              fg = CLIP_TUs[1:400], bg = bg5, tss_mode = "main", 
              width = 100, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4)



unlist(GRangesList(systemTUs[bidir$NDR_tag]$ohler))

viewTSS(systemTUs[div_sensitive$NDR_tag], width = 200, shift_val = -10, tss_mode = "main")

as.data.frame(div) %>%
  ggplot(aes(x = iscoding, y = logFC)) + 
  geom_boxplot(aes(color=factor(iscoding)), outlier.shape=NA, width = 0.5, lwd=1.25, notch = T) + 
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()

as.data.frame(div) %>%
  ggplot(aes(x = logFC, y = -log10(FDR))) + 
  geom_point(aes(color=factor(iscoding))) + 
  theme_minimal()

#mate-wise
check = as.data.frame(div) %>% 
  group_by(NDR) %>% 
  arrange(iscoding, .by_group = T) %>% 
  mutate(diff = abs(lag(logFC)/logFC))


##toy motif search
df_res = streme_wrapper(testreg1, testreg2, write_dir = "Shannahan/test")

#plotting the positional distribution
matplot(as.numeric(str_split(trimws(df_res$site_distr[13], "left"), " ")[[1]]), type = "l")








