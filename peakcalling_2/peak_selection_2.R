source("/Users/IEO5559/enhancedClip/Shannahan/feat_extr_kmer_functions.R")
source("/Users/IEO5559/enhancedClip/functions_2.R")
source("/Users/IEO5559/enhancedClip/Shannahan/functions.R")
library(seqinr)
library(LncFinder)
options(meme_bin = "/Users/IEO5559/bin/bin")
genome = BSgenome.Hsapiens.UCSC.hg38.masked

#slop
#top1000_slop100 = slop(clipper_peaks[1:1000], ext = 100)


streme_wrapper(fg_regions = slop(clipper_peaks, ext = 100), 
               bg_regions = bg3_init, order = 2, 
               write_dir = "streme_searches/vs_bg3/", maxw = 6)

streme_wrapper(fg_regions = slop(clipper_peaks, ext = 100), 
               bg_regions = clipper_peaks, order = 2, 
               write_dir = "streme_searches/slop100/", maxw = 6, nref = 20)


###### adjust streme params for RBP motif
streme_wrapper(fg_regions = slop(clipper_peaks, ext = 100), 
               bg_regions = clipper_peaks, order = 1, 
               write_dir = "streme_searches/slop100_subsearches/asdfdsf/", 
               minw = 4, maxw = 6, seed = 111, 
               nref = 5, neval = 25, holdout = 0.2, niter = 20)

streme_wrapper(fg_regions = slop(clipper_peaks, ext = 100), 
               bg_regions = clipper_peaks, order = 2, 
               write_dir = "streme_searches/slop100_subsearches/asdfdsf2/", 
               minw = 4, maxw = 6, seed = 111, 
               nref = 10, neval = 25, holdout = 0.2, niter = 1)

streme_wrapper(fg_regions = slop(clipper_peaks, ext = 100), 
               bg_regions = clipper_peaks, order = 2, 
               write_dir = "streme_searches/slop100_subsearches/asdfdsf3/", 
               minw = 4, maxw = 6, seed = 1511, 
               nref = 10, neval = 25, holdout = 0.2, niter = 1)

streme_wrapper(fg_regions = slop(clipper_peaks, ext = 100), 
               bg_regions = clipper_peaks, order = 2, 
               write_dir = "streme_searches/slop100_subsearches/asdfdsf4/", 
               minw = 4, maxw = 6, seed = 1511, 
               nref = 10, neval = 1, holdout = 0.2, niter = 1)








