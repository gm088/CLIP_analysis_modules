source("/Users/IEO5559/enhancedClip/Shannahan/feat_extr_kmer_functions.R")
source("~/enhancedClip/functions_2.R")
source("/Users/IEO5559/enhancedClip/Shannahan/functions.R")
options(meme_bin = "/Users/IEO5559/meme/bin")

library(seqinr)
library(LncFinder)

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

extension = 1000

TSS_and_motif(outdir = "outputs/bg2", filename = "CLIP_TU_vs_bg2.pdf", 
              fg = CLIP_TUs, bg = bg2, tss_mode = "main", 
              width = extension, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4,
              fg_name = "CLIP_TUs", bg_name = "coding unaffected")

TSS_and_motif(outdir = "outputs/bg3", filename = "CLIP_TU_vs_bg3.pdf", 
              fg = CLIP_TUs, bg = bg3, tss_mode = "main", 
              width = extension, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4,
              fg_name = "CLIP_TUs", bg_name = "divergent units (unaffected)")

TSS_and_motif(outdir = "outputs/tophalf_vs_bg3", filename = "tophalf_vs_bg3.pdf", 
              fg = CLIP_TUs[1:400], bg = bg3, tss_mode = "main", 
              width = extension, shift_val = -2, minw = 4, maxw = 20, streme_order = 1:4,
              fg_name = "CLIP_TUs top half", bg_name = "divergent units (unaffected)")

TSS_and_motif(outdir = "outputs/bg5", filename = "CLIP_TU_vs_bg5.pdf", 
              fg = CLIP_TUs, bg = bg5, tss_mode = "main", 
              width = extension, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4,
              fg_name = "CLIP_TUs", bg_name = "mates of known pa-RNA")

TSS_and_motif(outdir = "outputs/bg5_pro", filename = "CLIP_TU_vs_bg5.pdf", 
              fg = CLIP_TUs, bg = bg5, tss_mode = "main", usePROSEQ = T,
              width = extension, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4,
              fg_name = "CLIP_TUs", bg_name = "mates of sensitive pa-RNA")

TSS_and_motif(outdir = "outputs/tophalf_vs_bg3_pro", filename = "tophalf_vs_bg3_pro.pdf", 
              fg = CLIP_TUs[1:400], bg = bg3, tss_mode = "main", usePROSEQ = T,
              width = extension, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4,
              fg_name = "CLIP_TUs top half", bg_name = "divergent units (unaffected)")

TSS_and_motif(outdir = "outputs/div_sensitive", filename = "div_vs_div_sensitive.pdf", 
              fg = systemTUs[div_sensitive$NDR_tag], bg = bg5, tss_mode = "main", 
              width = extension, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4,
              fg_name = "sensitive pa-RNA", bg_name = "mates of sensitive pa-RNA")

TSS_and_motif(outdir = "outputs/tophalf_vs_bg5", filename = "tophalf_vs_bg5.pdf", 
              fg = CLIP_TUs[1:400], bg = bg5, tss_mode = "main", 
              width = extension, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4,
              fg_name = "CLIP_TUs top half", bg_name = "mates of sensitive pa-RNA")


#####testing
TSS_and_motif(outdir = "outputs/vs_bg5_rightmost", filename = "vs_bg5_rightmost.pdf", 
              fg = CLIP_TUs, bg = bg5, fg_tss_mode = "main", bg_tss_mode = "rightmost",
              width = 100, shift_val = -2, minw = 4, maxw = 30, streme_order = 1:4,
              fg_name = "CLIP_TUs top half", bg_name = "mates of sensitive pa-RNA")

TSS_and_motif(outdir = "outputs/div_sensitive_rightmost", filename = "div_vs_div_sensitive.pdf", 
              fg = systemTUs[div_sensitive$NDR_tag], bg = bg5, fg_tss_mode = "rightmost", bg_tss_mode = "rightmost",
              width = 200, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4,
              fg_name = "sensitive pa-RNA", bg_name = "mates of sensitive pa-RNA")

TSS_and_motif(outdir = "outputs/bidir_left_vs_right", filename = "bidir_left_vs_right", 
              fg = systemTUs[bidir$NDR_tag], bg = systemTUs[bidir$NDR_tag], fg_tss_mode = "leftmost", bg_tss_mode = "rightmost",
              width = 200, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4,
              fg_name = "bidir leftmost", bg_name = "bidir rightmost")

TSS_and_motif(outdir = "outputs/tophalf_vs_bg5_meme", filename = "tophalf_vs_bg5.pdf", 
              fg = CLIP_TUs[1:400], bg = bg5,  fg_tss_mode = "main", bg_tss_mode = "main",
              width = 100, shift_val = -2, minw = 4, maxw = 10, streme_order = 1:4,
              fg_name = "CLIP_TUs top half", bg_name = "mates of sensitive pa-RNA", usestreme = F)


###### meme test

fg = getSeq(genome, shiftStranded(resize(get_uniq_TSS(CLIP_TUs, mode = "main"), width = 100, fix = "start"), -2))
names(fg) = get_uniq_TSS(CLIP_TUs, mode = "main")$tag
bg = getSeq(genome, shiftStranded(resize(get_uniq_TSS(bg5, mode = "main"), width = 100, fix = "start"), -2))
names(bg) = get_uniq_TSS(bg5, mode = "main")$tag

results = runMeme(
  fg,
  control = bg,
  outdir = "Shannahan/outputs/tophalf_vs_bg5_meme",
  alph = "rna",
  objfun = "de",
  minw = 4,
  maxw = 10,
  nmotifs = 30,
  combined_sites = T,
  shuf = 2,
  parse_genomic_coord = F,
  mod = "anr"
)


####RNA fold test 
width=20
writeXStringSet(
  getSeq(genome, shiftStranded(resize(get_uniq_TSS(CLIP_TUs, mode = "main"), width = width, fix = "start"), -2)), 
  filepath = "/Users/ieo5559/tmpfastas/fg.fa")
writeXStringSet(
  getSeq(genome, shiftStranded(resize(get_uniq_TSS(bg5, mode = "main"), width = width, fix = "start"), -2)), 
  filepath = "/Users/ieo5559/tmpfastas/bg.fa")

fg_seqs = read.fasta("/Users/ieo5559/tmpfastas/fg.fa")
bg_seqs = read.fasta("/Users/ieo5559/tmpfastas/bg.fa")

folded_fg = run_RNAfold(fg_seqs)
folded_bg = run_RNAfold(bg_seqs)

hist(as.numeric(unlist(lapply(folded_bg, function(x){x[[3]]}))))
hist(as.numeric(unlist(lapply(folded_fg, function(x){x[[3]]}))), add = T)

lapply(folded_fg, function(x){x[2]})[1:5]

system2("RNAfold", args = c("-i /Users/ieo5559/tmpfastas/fg.fa", "-p"))


####write fasta files for other motif search tools

ups = 20
downs = 150
dir1 = "/Users/ieo5559/enhancedClip/Shannahan/region_fastas/5primeGRO/"
dir2 = "/Users/ieo5559/enhancedClip/Shannahan/region_fastas/PROSEQ/"

writeregions(systemtus_list = list(CLIP_TUs, bg3, bg5), ups = ups, downs = downs, 
             tags = list("CLIP_TU", "bg3", "bg5"), dir = dir1, PROSEQ = F)

writeregions(systemtus_list = list(CLIP_TUs, bg3, bg5), ups = ups, downs = downs, 
             tags = list("CLIP_TU_pro", "bg3_pro", "bg5_pro"), dir = dir2, PROSEQ = T)

####

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








