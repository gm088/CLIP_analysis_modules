library(csaw)
library(rtracklayer)
library(edgeR)
library(BiocParallel)
library(GenomicRanges)
source("utility_functions.R")

##### INITIALISE ##### 

bams_f = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_fwd_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/94_SE_fwd_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/95_SE_fwd_clean.bam")

bams_r = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_rev_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/94_SE_rev_clean.bam",
           "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/95_SE_rev_clean.bam")

samples = 1:3
Desc = c("INPUT", "IP2", "IP2")
metadata = data.frame(samples, Desc)
desmat <- model.matrix(~factor(metadata$Desc, levels = c("INPUT", "IP2")))
colnames(desmat) = c("Intercept", "IP2")

blacklist = import("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/blacklists/my_new_clip_blacklist_1539_6Dec_2022.bed")
win_width = 250
min_reads_per_win = 50
win_spacing = 250
param = readParam(minq=10, pe="none", discard = blacklist) # important - reads extracted from bams with these params

##### COUNTING ##### 

demo_plus <- windowCounts(bams_f, param=param, BPPARAM=MulticoreParam(), width = win_width, filter=min_reads_per_win, spacing = win_spacing, ext = NA)

demo_rev <- windowCounts(bams_r, param=param, BPPARAM=MulticoreParam(), width = win_width, filter=min_reads_per_win, spacing = win_spacing, ext = NA)

strand(rowRanges(demo_plus)) = "+"
strand(rowRanges(demo_rev)) = "-"

##### SIZE FACTORS ##### 

binned_f <- windowCounts(bams_f, bin=TRUE, width=10000, param=param, BPPARAM=MulticoreParam())
binned_r <- windowCounts(bams_r, bin=TRUE, width=10000, param=param, BPPARAM=MulticoreParam())

demo_plus = normFactors(binned_f, se.out=demo_plus, assay.id=1)
demo_rev = normFactors(binned_r, se.out=demo_rev, assay.id=1)

#import restr regions

#restr = import("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mysiW/mysiW_slopped_200bp.bed")

#demo_plus_restr = demo_plus[rowRanges(demo_plus) %in% subsetByOverlaps(rowRanges(demo_plus), restr, ignore.strand = F)]
#demo_rev_restr = demo_rev[rowRanges(demo_rev) %in% subsetByOverlaps(rowRanges(demo_rev), restr, ignore.strand = F)]

#In edgeR, the log-transformed overall NB mean is referred to as the average abundance.
#This is computed with the aveLogCPM() function, as shown below for each window

#### FITTING

y_f = asDGEList(demo_plus, assay.id=1)
y_f = estimateDisp(y_f,desmat)
summary(y_f$trended.dispersion)
fit_f = glmQLFit(y_f, desmat, robust=T)
summary(fit_f$var.post)

#pdf("dispersions_plus.pdf")
par(mfrow=c(1,2))
o <- order(y_f$AveLogCPM)
plot(y_f$AveLogCPM[o], sqrt(y_f$trended.dispersion[o]), type="l", lwd=2,
     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
     ylab=("Biological coefficient of variation"))
plotQLDisp(fit_f)
#dev.off()

results_f = glmQLFTest(fit_f, coef="IP2")
#results_f = glmQLFTest(fit_f, contrast = c(0,-1,1))

y_r = asDGEList(demo_rev, assay.id=1)
y_r = estimateDisp(y_r,desmat)
summary(y_f$trended.dispersion)
fit_r = glmQLFit(y_r, desmat, robust=T)
summary(fit_r$var.post)

#pdf("dispersions_minus.pdf")
par(mfrow=c(1,2))
o <- order(y_r$AveLogCPM)
plot(y_r$AveLogCPM[o], sqrt(y_r$trended.dispersion[o]), type="l", lwd=2,
     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
     ylab=("Biological coefficient of variation"))
plotQLDisp(fit_r)
#dev.off()

results_r = glmQLFTest(fit_r, coef="IP2")
#results_r = glmQLFTest(fit_r, contrast = c(0,-1,1))

#### APPENDING INFO

mcols(rowRanges(demo_plus)) = cbind(mcols(rowRanges(demo_plus)), results_f$table)
mcols(rowRanges(demo_rev)) = cbind(mcols(rowRanges(demo_rev)), results_r$table)

####  choose Pval 0.001, log2FC > 2
pvalthresh = 1e-2
fcthresh = 2

goodwins = unlist(GRangesList(rowRanges(demo_plus)[rowRanges(demo_plus)$logFC > fcthresh & rowRanges(demo_plus)$PValue < pvalthresh], 
                              rowRanges(demo_rev)[rowRanges(demo_rev)$logFC > fcthresh & rowRanges(demo_rev)$PValue < pvalthresh]))

#allwins = unlist(GRangesList(rowRanges(demo_plus), rowRanges(demo_rev)))
#mcols(allwins) = cbind(mcols(allwins), rbind(assay(demo_plus), assay(demo_rev)))

#merge - do not set a max width or whatever, I don't want overlapping peaks

merged = mergeWindows(goodwins, max.width = NULL, tol=2000, ignore.strand = F) #2nd is the max dist between adjacent windows
merged$regions

summary(width(merged$regions))

tabcom <- combineTests(merged$ids, mcols(goodwins))
mcols(merged$regions) = tabcom
is.sig.region = tabcom$FDR <= 1e-5
summary(is.sig.region)


#check
merged$regions[merged$regions$FDR <= 1e-5 & merged$regions$direction=="up"]

fdrsig = merged$regions[merged$regions$FDR <= 1e-5 & merged$regions$direction=="up"]

df <- data.frame(seqnames=seqnames(fdrsig),
                 starts=start(fdrsig),
                 ends=end(fdrsig),
                 scores=fdrsig$FDR,
                 FC=fdrsig$rep.logFC,
                 strand=strand(fdrsig))

write.table(df, "with_clean_bams_12Dec_2022_no_tRNA_final_win150bp_spacing50bp.bed", quote = F, col.names = F, row.names = F, sep = "\t")

#save.image(file = "simple_desmat_w_empty_workspace_postblacklist.RData")

#2Dec 2022
#save.image(file = "with_clean_bams_1742_2Dec_2022.RData")

#5Dec 2022
save.image(file = "with_clean_bams_12Dec_2022_no_tRNA_final_win150bp_spacing50bp.RData")

######################## 29th Nov 2022

#savetoblacklist = merged$regions[tabcom$FDR < 0.0001,]

#bl <- data.frame(seqnames=seqnames(savetoblacklist),
#                 starts=start(savetoblacklist),
#                 ends=end(savetoblacklist),
#                 scores=savetoblacklist$FDR,
#                 FC=savetoblacklist$rep.logFC,
#                 strand=strand(savetoblacklist))
#
#write.table(bl, "my_clip_blacklist.bed", quote = F, col.names = F, row.names = F, sep = "\t")

####################### 30th Nov 2022

# savetoblacklist = unlist(GRangesList(blacklist, merged$regions[merged$regions$FDR < 0.0009]))
# 
# bl <- data.frame(seqnames=seqnames(savetoblacklist),
#                  starts=start(savetoblacklist),
#                  ends=end(savetoblacklist))
# 
# write.table(bl, "my_clip_blacklist.bed", quote = F, col.names = F, row.names = F, sep = "\t")
# 
#results_all = rbind(cbind(results_f$table, "deviance" = results_f$deviance, "dispersion" = results_f$dispersion, "var.prior" = results_f$var.prior, "residual" = results_f$df.residual), 
#                    cbind(results_r$table, "deviance" = results_r$deviance, "dispersion" = results_r$dispersion, "var.prior" = results_r$var.prior, "residual" = results_r$df.residual))
#results_notM = results_all[as.vector(seqnames(allwins) != "chrM"),]
#notM = allwins[seqnames(allwins) != "chrM"]
#
#mcols(notM) = cbind(mcols(notM), results_notM)

####################### 1st Dec 2022

# surrounds = 1000
# sus = merged$regions[merged$regions$FDR < 0.001]
# neigh = suppressWarnings(resize(sus, surrounds, fix="center"))
# 
# neigh_f = neigh[strand(neigh)=="+"]
# sus_f = sus[strand(sus)=="+"]
# goodwins_f = goodwins[strand(goodwins)=="+"]
# wider_f = regionCounts(bams_f, neigh_f , param = param)
# not_wider_f = regionCounts(bams_f, sus_f, param = param)
# 
# neigh_r = neigh[strand(neigh)=="-"]
# sus_r = sus[strand(sus)=="-"]
# goodwins_r = goodwins[strand(goodwins)=="-"]
# wider_r = regionCounts(bams_r, neigh_r , param = param)
# not_wider_r = regionCounts(bams_r, sus_r, param = param)
# 
# #### find those regions for which neighbourhood has zero reads (i.e. abs(not_wider - wider) == 0 (or v small))
# 
# diff_f = abs(assay(not_wider_f) - assay(wider_f))
# inds_f = rowSums(diff_f) == 0
# tower_p = sus_f[inds_f]
# 
# diff_r = abs(assay(not_wider_r) - assay(wider_r))
# inds_r = rowSums(diff_r) == 0
# tower_r = sus_f[inds_r]
# 
# almost_newtowers = unlist(GRangesList(tower_p, tower_r))
# newtowers = almost_newtowers[almost_newtowers$num.tests < 3]
# addtobl = newtowers[newtowers$FDR < 0.00001]
# 
# newbl = reduce(unlist(GRangesList(addtobl, blacklist)))
# export(newbl, "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/my_new_clip_blacklist_1742_2Dec_2022.bed")


########6th Dec 2022

#count to distinguish towers - final, I hope
# question = fdrsig[fdrsig$FDR < 1e-7]
# question = resize(question, width = 50, fix = "center")
# q_neigh = resize(question, width = 150, fix = "center")
# 
# question_plus = question[strand(question) == "+"]
# question_minus = question[strand(question) == "-"]
# q_neigh_plus = q_neigh[strand(q_neigh) == "+"]
# q_neigh_minus = q_neigh[strand(q_neigh) == "-"]
# 
# counts_plus = regionCounts(bams_f, question_plus , param = param)
# counts_minus = regionCounts(bams_r, question_minus , param = param)
# counts_neigh_plus = regionCounts(bams_f, q_neigh_plus , param = param)
# counts_neigh_minus = regionCounts(bams_r, q_neigh_minus , param = param)
# 
# diff_plus = assay(counts_neigh_plus)[,c(3,4)] - assay(counts_plus)[,c(3,4)]
# diff_minus = assay(counts_neigh_minus)[,c(3,4)] - assay(counts_minus)[,c(3,4)]
# diff = rbind(diff_plus, diff_minus)
# 
# mcols(question) = cbind(mcols(question), diff)
# 
# towers = question[rowSums(diff) < 5]

#check = merged$regions[merged$regions$FDR < 1e-8]
#heck2 = merged$regions[merged$regions$FDR < 1e-6]

#########polishing - 22 Dec 2022

source("utility_functions.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")
library(GenomicAlignments)

sam_f = c("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/IP7_DCIP_SE_fwd.bam")
sam_r = c("/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/IP7_DCIP_SE_rev.bam")

#check the profiles
tf = merged$regions[width(merged$regions) == 250]
names(tf) = paste("id", 1:length(tf), sep = "")
plus_prof = coverageFromFile(txdb = NULL, file = sam_f, gr = tf[strand(tf)=="+"], change_chr_names = F, strand = "+")
minus_prof = coverageFromFile(txdb = NULL, file = sam_r, gr = tf[strand(tf)=="-"], change_chr_names = F, strand = "-")

profiles = rbind(do.call('rbind', plus_prof), do.call('rbind', minus_prof))
norm = nainorm(profiles)

##gradient - each profile is 250bp

ddnorm = lapply(1:length(tf), function(x){norm[x,2:250] - norm[x, 1:249]})

#given that the max (of normalised profiles) is 1, profiles in which max of gradient is close to 1 is a tower

maxgrads = unlist(lapply(ddnorm, max))
bad_indices = maxgrads > 0.75

check = norm[bad_indices,]

##just add these regions to the blacklist and call it a day






