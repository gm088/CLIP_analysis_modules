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

#I want the procap, IP7, DCIP, SMInput, PROSEQ
b_f = c("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/14_r2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/94_merged_read2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/95_r2_for.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/LisCap/bigwigs/SRR22522159_fwd_5p.bw")

b_r = c("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/14_r2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/cleanbams/94_merged_read2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2bw/95_r2_rev.bw",
        "/hpcnfs/data/GN2/gmandana/take2project/LisCap/bigwigs/SRR22522159_rev_5p.bw")


# ### we need to find all TSS custers using the procap bigwigs
# 
# procap_plus = import("procap/procap_for.bw", as = 'GRanges'); strand(procap_plus) = "+"
# procap_minus = import("procap/procap_rev.bw", as = 'GRanges'); strand(procap_minus) = "-"
# procap_minus$score = abs(procap_minus$score)
# 
# # choose a threshold
# thresh = 10
# 
# alltss = c(procap_plus[procap_plus$score > thresh], 
#            procap_minus[procap_minus$score > thresh])
# 
# ###clustering
# clus = csaw::mergeWindows(alltss, ignore.strand = F,
#                           tol = 20)
# 
# alltss$clus = clus$ids
# 
# clus_scores = as.data.frame(alltss) %>% 
#   group_by(clus) %>% 
#   dplyr::select(c("score", "clus")) %>% 
#   summarise(clus_score = sum(score))
# 
# clus$regions$score = clus_scores$clus_score
# 
# clustered_tss = clus$regions


################## 19Dec - from macs clsuters to point source

### for each cluster, just get the max.
macspeaks = customimportbed6("/hpcnfs/data/GN2/gmandana/take2project/LisCap/macs/round1/MACS_peaks.narrowPeak")
macspeaks$id = paste0("id", 1:length(macspeaks))
bam_f = "/hpcnfs/data/GN2/gmandana/take2project/LisCap/mapped/SRR22522159/fwd_5prime_sorted.bam"
bam_r = "/hpcnfs/data/GN2/gmandana/take2project/LisCap/mapped/SRR22522159/rev_5prime_sorted.bam"

mpp = macspeaks[strand(macspeaks) == "+"]
mpm = macspeaks[strand(macspeaks) == "-"]

cov_f = coverage(bam_f)
cov_r = coverage(bam_r)
#which.max(as.vector(cov_f[macspeaks[strand(macspeaks) == "+"]][[1]]))

offset_p = sapply(cov_f[mpp], function(x){which.max(as.vector(x))})
offset_m = sapply(cov_r[mpm], function(x){which.max(as.vector(x))})

#make new granges

new = c(GRanges(seqnames = seqnames(mpp), ranges = IRanges(start = start(mpp)+offset_p), strand = "+"),
GRanges(seqnames = seqnames(mpm), ranges = IRanges(start = end(mpm)-(width(mpm)-offset_m), width = 1), strand = "-"))

export(new, "PROCAP_maxima_rep1_v1.bed")


########### 29 Jan 2024 - define procap regions for metaplot

### peak cluster level
peak_clusters = readRDS("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_26Oct/clusters_final.RDS")
clus_ranges = makeGRangesFromDataFrame(peak_clusters, keep.extra.columns = T)
#procap_point = new
procap_point = import("PROCAP_maxima_rep1_v1.bed")
dis = as.data.frame(
  distanceToNearest(promoters(clus_ranges, upstream = 1, downstream = 0), procap_point, ignore.strand = F)
)

maxdist = 1000

clus_ranges$disttoprocap = dis$distance
goodclusranges = clus_ranges[dis[dis$distance < maxdist,]$queryHits] # good ones

#these are your ranges to plot at cluster level
goodclusranges_procap = procap_point[dis[dis$distance < maxdist,]$subjectHits]

procap_aligned_prof = getcovmat(list(b_f, b_r), 
                                promoters(goodclusranges_procap, upstream = 100, downstream = 1000), 
                                nbins = 110, norm = T)

matplot(t(procap_aligned_prof), type = "l")


### siW...
siwregions = mysiW[grepl(mysiW$name, pattern = "pa-RNA")]
dis2siw = as.data.frame(
  distanceToNearest(promoters(siwregions, upstream = 1, downstream = 0), procap_point, ignore.strand = F)
)

goodsiwranges_procap = procap_point[dis2siw[dis2siw$distance < maxdist,]$subjectHits]

procap_aligned_prof_siw = getcovmat(list(b_f, b_r), 
                                promoters(goodsiwranges_procap, upstream = 100, downstream = 1000), 
                                nbins = 110, norm = F)

matplot(t(procap_aligned_prof_siw), type = "l")



#### multiple peaks per transcript
clipper_peaks[duplicated(clipper_peaks$V4)]





########## procap - use V10 as offset (see narrowpeak formatpctx)

#make new granges

mpp = macspeaks[strand(macspeaks) == "+"]
mpm = macspeaks[strand(macspeaks) == "-"]

offset_p = mpp$V10
offset_m = mpm$V10

new2p = GRanges(seqnames = seqnames(mpp), ranges = IRanges(start = start(mpp)+offset_p), strand = "+")
new2m = GRanges(seqnames = seqnames(mpm), ranges = IRanges(start = end(mpm)-(width(mpm)-offset_m), width = 1), strand = "-")

new2p$score = mpp$V9
new2m$score = mpm$V9
new2 = c(new2p, new2m)

##append scores

export(new2, "PROCAP_maxima_rep1_v2.bed")



















