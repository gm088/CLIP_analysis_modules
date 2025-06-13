source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")

####given a set of foreground sequences, what distinguishes it from a background?

# set pseudorandom number generator
set.seed(10)
options(meme_bin = "/hpcnfs/home/ieo5559/meme/bin")

genome = BSgenome.Hsapiens.UCSC.hg38.masked
peaks_file = "new_top_peaks_clipper.bed"
bgfile = "all5prime_expressed_unique_500bp_notdiffexpr_notbound.bed"
bgfile2 = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/files_for_svm/proseq_position_calibrated.bed"
motifs = readRDS("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/files_for_svm/merged_motifs_2.RDS")
repeats = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/files_for_svm/repeats_families.bed")
simple_repeats = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/files_for_svm/simple_repeats_fixed.bed")
fimo_bg_file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/files_for_svm/bg_500ups1kbdowns.bg"

###train

fg = customimportbed6(peaks_file)
bg = import(bgfile)
bg_pro = customimportbed6(bgfile2)

inputs = collect(fg, bg, motifs, repeats, simple_repeats, expand = F,
                 500, 1000, fimo_bg_file, 1e-3, 1000000, list(h3k79me2, "h3k79me2"), 
                 list(h3k27ac, "h3k27ac"))
inputs_pro = collect(fg, bg_pro, motifs, repeats, simple_repeats, expand = F,
                     500, 1000, fimo_bg_file, 1e-3, 100000, list(h3k79me2, "h3k79me2"), 
                     list(h3k27ac, "h3k27ac"))

saveRDS(inputs, file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/files_for_svm/inputs_6Jun_346mer.RDS")
saveRDS(inputs_pro, file = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/files_for_svm/inputs_pro_6Jun_346mer.RDS")
