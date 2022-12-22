library(csaw)
library(rtracklayer)
library(edgeR)
library(BiocParallel)
library(GenomicRanges)
source("utility_functions.R")

##### INITIALISE ##### 

bams_f = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/94_SE_fwd.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/95_SE_fwd.bam")

bams_r = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/94_SE_rev.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/95_SE_rev.bam")

samples = 1:3
Desc = c("INPUT", "IP2", "IP2")
metadata = data.frame(samples, Desc)
desmat <- model.matrix(~factor(metadata$Desc, levels = c("INPUT", "IP2")))
colnames(desmat) = c("Intercept", "IP2")

blacklist = import("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/my_clip_blacklist.bed")
win_width = 50
min_reads_per_win = 20
win_spacing = 25
param = readParam(minq=10, pe="none", discard = blacklist) # important - reads extracted from bams with these params

##### COUNTING ##### 

demo_plus <- windowCounts(bams_f, param=param, BPPARAM=MulticoreParam(), width = win_width, filter=min_reads_per_win, spacing = win_spacing, ext = NA)
demo_rev <- windowCounts(bams_r, param=param, BPPARAM=MulticoreParam(), width = win_width, filter=min_reads_per_win, spacing = win_spacing, ext = NA)

strand(rowRanges(demo_plus)) = "+"
strand(rowRanges(demo_rev)) = "-"

##### neighbourhood count for towers ##### 

surrounds <- 1000
neighbor_f <- suppressWarnings(resize(rowRanges(demo_plus), surrounds, fix="center"))
wider_f <- regionCounts(bams_f, regions=neighbor_f, param=param)
filter.stat_f <- filterWindowsLocal(demo_plus, wider_f, assay.data = 1, assay.back = 1)

summary(filter.stat_f$filter)
keep_f <- filter.stat_f$filter > log2(3)
sum(keep_f)

hist(filter.stat_f$filter, xlab="Log-fold change from local background", 
     breaks=100, main="", col="grey80", xlim=c(0, 5))
abline(v=log2(3), col="red", lwd=2)

###minus

neighbor_r <- suppressWarnings(resize(rowRanges(demo_rev), surrounds, fix="center"))
wider_r <- regionCounts(bams_r, regions=neighbor_r, param=param)
filter.stat_r <- filterWindowsLocal(demo_rev, wider_r)

summary(filter.stat_r$filter)
keep_r <- filter.stat_r$filter > log2(3)
sum(keep_r)

hist(filter.stat_r$filter, xlab="Log-fold change from local background", 
     breaks=100, main="", col="grey80", xlim=c(0, 5))
abline(v=log2(3), col="red", lwd=2)

save.image(file = "finding_clip_towers.RData")

c1 = demo_plus[rowSums(assay(demo_plus)) > 100]
c2 = wider_f[rowSums(assay(demo_plus)) > 100]

check = (assay(c1)+1)/(assay(c2)+1)






