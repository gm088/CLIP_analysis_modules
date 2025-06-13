library(csaw)
library(rtracklayer)
library(edgeR)
library(BiocParallel)
library(GenomicRanges)

initialise = function(){
  bams_f <<- c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/71_SE_fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/73_SE_fwd.bam")
  
  bams_r <<- c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/71_SE_rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/73_SE_rev.bam")
  
  samples = 1:3
  Desc = c("INPUT","ZC3H4","ZC3H4")
  metadata = data.frame(samples, Desc)
  desmat <<- model.matrix(~factor(metadata$Desc))
  colnames(desmat) = c("Intercept", "ZC3H4")
  desmat <<- desmat
}

windowcount = function(bf, br, param, win_width, min_reads_per_win, win_spacing){
  
  demo_plus <<- windowCounts(bf, param=param, BPPARAM=MulticoreParam(), width = win_width, filter=min_reads_per_win, spacing = win_spacing, ext = NA)
  print("done forward")
  demo_rev <<- windowCounts(br, param=param, BPPARAM=MulticoreParam(), width = win_width, filter=min_reads_per_win, spacing = win_spacing, ext = NA)
  print("done reverse")
  strand(rowRanges(demo_plus)) = "+"
  strand(rowRanges(demo_rev)) = "-"
  demo_plus <<- demo_plus
  demo_rev <<- demo_rev
  
}

getAbundances = function(rsexp){
  
  abundances <- aveLogCPM(asDGEList(rsexp, assay.id=1))
  print(summary(abundances))
  #return(abundances)
  
}


visdisp = function(y, fit){
  
  par(mfrow=c(1,2))
  o <- order(y$AveLogCPM)
  plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
       ylim=c(0, 5), xlab=expression("Ave."~Log[2]~"CPM"),
       ylab=("Biological coefficient of variation"))
  plotQLDisp(fit)
  
}























