library(csaw)
library(rtracklayer)
library(edgeR)
library(BiocParallel)
library(GenomicRanges)
source("utility_functions.R")

##### INITIALISE ##### 

initialise()

win_width = 50
param = readParam(minq=10, pe="none") # important - reads extracted from bams with these params

##### COUNTING ##### 

windowcount(bf = bams_f, br = bams_r, param = param, win_width = win_width, min_reads_per_win = 5, win_spacing = 25)

##### SIZE FACTORS ##### 

binned_f <- windowCounts(bams_f, bin=TRUE, width=20000, param=param, BPPARAM=MulticoreParam())
binned_r <- windowCounts(bams_r, bin=TRUE, width=20000, param=param, BPPARAM=MulticoreParam())

demo_plus = normFactors(binned_f, se.out=demo_plus, assay.id=1)
demo_rev = normFactors(binned_r, se.out=demo_rev, assay.id=1)

#### count-based filtering. ####

keep_f = rowSums(assay(demo_plus)) >= 20   # at least 20 reads across all samples in a window
keep_r = rowSums(assay(demo_rev)) >= 20

#import restr regions

#restr = import("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mysiW/mysiW_slopped_200bp.bed")

#demo_plus_restr = demo_plus[rowRanges(demo_plus) %in% subsetByOverlaps(rowRanges(demo_plus), restr, ignore.strand = F)]
#demo_rev_restr = demo_rev[rowRanges(demo_rev) %in% subsetByOverlaps(rowRanges(demo_rev), restr, ignore.strand = F)]

#### it's ok for norm factors, no composition bias

#In edgeR, the log-transformed overall NB mean is referred to as the average abundance.
#This is computed with the aveLogCPM() function, as shown below for each window.

abundances_demo_plus = getAbundances(demo_plus[keep_f])
abundances_demo_rev = getAbundances(demo_rev[keep_r])

demo_plus_filt = demo_plus[keep_f]
demo_rev_filt = demo_rev[keep_r]

#demo_plus_filt = normFactors(demo_plus_filt, assay.id = 1, se.out = demo_plus_filt)
#demo_rev_filt = normFactors(demo_rev_filt, assay.id = 1, se.out = demo_rev_filt)

###new lib sizes?

regcounts_plus = regionCounts(bams_f, reduce(rowRanges(demo_plus_filt)), param=param)
regcounts_rev = regionCounts(bams_r, reduce(rowRanges(demo_rev_filt)), param=param)

newlibsizes_plus = colSums(assay(regcounts_plus))
newlibsizes_rev = colSums(assay(regcounts_rev))

demo_plus_filt$totals = newlibsizes_plus
demo_rev_filt$totals = newlibsizes_rev

rm(regcounts_plus, regcounts_rev, binned_f, binned_r)

### new norm factors?

demo_plus_filt = normFactors(demo_plus_filt, assay.id = 1, se.out = demo_plus_filt)
demo_rev_filt = normFactors(demo_rev_filt, assay.id = 1, se.out = demo_rev_filt)

#### FITTING

y_f = asDGEList(demo_plus_filt, assay.id=1)
y_f = estimateDisp(y_f,desmat)
summary(y_f$trended.dispersion)
fit_f = glmQLFit(y_f, desmat, robust=T)
summary(fit_f$var.post)

visdisp(y_f, fit_f)
results_f = glmQLFTest(fit_f, coef="ZC3H4")

y_r = asDGEList(demo_rev_filt, assay.id=1)
y_r = estimateDisp(y_r,desmat)
summary(y_f$trended.dispersion)
fit_r = glmQLFit(y_r, desmat, robust=T)
summary(fit_r$var.post)

results_r = glmQLFTest(fit_r, coef="ZC3H4")

#### APPENDING INFO

mcols(rowRanges(demo_plus_filt)) = cbind(mcols(rowRanges(demo_plus_filt)), results_f$table)
mcols(rowRanges(demo_rev_filt)) = cbind(mcols(rowRanges(demo_rev_filt)), results_r$table)


####  choose Pval 0.001, log2FC > 2
pvalthresh = 0.005
fcthresh = 2

goodwins = unlist(GRangesList(rowRanges(demo_plus_restr)[rowRanges(demo_plus_restr)$logFC > fcthresh & rowRanges(demo_plus_restr)$PValue < pvalthresh], 
                              rowRanges(demo_rev_restr)[rowRanges(demo_rev_restr)$logFC > fcthresh & rowRanges(demo_rev_restr)$PValue < pvalthresh]))

allwins = unlist(GRangesList(rowRanges(demo_plus_restr), rowRanges(demo_rev_restr)))

#merge - do not set a max width or whatever, I don't want overlapping peaks

merged = mergeWindows(goodwins, max.width = NULL, tol=0, ignore.strand = F) #2nd is the max dist between adjacent windows
merged$regions

summary(width(merged$regions))

tabcom <- combineTests(merged$ids, mcols(goodwins))
is.sig.region <- tabcom$FDR <= 0.01
summary(is.sig.region)

mcols(merged$regions) = tabcom

#check
merged$regions[merged$regions$FDR <= 0.01 & merged$regions$direction=="up"]

fdrsig = merged$regions[merged$regions$FDR <= 0.01 & merged$regions$direction=="up"]




df <- data.frame(seqnames=seqnames(fdrsig),
                 starts=start(fdrsig),
                 ends=end(fdrsig),
                 scores=fdrsig$FDR,
                 FC=fdrsig$rep.logFC,
                 strand=strand(fdrsig))

write.table(df, "test_with_shift_win50_spacing25.bed", quote = F, col.names = F, row.names = F, sep = "\t")


