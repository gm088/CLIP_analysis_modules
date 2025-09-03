library(fitdistrplus)
library(BiocParallel)
library(forcats)
library(tools)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(abind)

function_to_match_expr = function(fg, bg, grain = 0.1, ssize = 2000){
  
  print("make sure fg rpkm in col 'V7' and bg rpkm in col V5")
  pc = 0.001
  
  #we will fit the rpkm of the foreground seqs to a lognormal dist
  #dis$estimate[[1]] is the mean and [[2]] is the sd
  dis = fitdist(fg$V7+pc, distr = "lnorm", method = "mme")
  lmean = dis$estimate[[1]]
  lsd = dis$estimate[[2]]
  #pseudocount
  #construct probability vector for sampling
  #use the desnisty function to estimate the probability of each rpkm value
  #each rpkm value in a window of 'grain'
  
  probs = plnorm(bg_coding$V5+pc+grain, meanlog = lmean, sdlog = lsd) - 
    plnorm(bg_coding$V5+pc, meanlog = lmean, sdlog = lsd)
  
  #hist(bg_coding[sample(1:length(bg_coding), size = 2000, replace = F, prob = probs)]$V5+pc, 
  #     breaks = 100,
  #     col=rgb(1,0,0,1/4))
  #hist(fg$V7+pc, breaks = 100, add = T, col=rgb(0,0,1,1/4))
  #hist(bg_coding$V5+pc, breaks = 10000, add = T, col=rgb(0,1,0,1/4))
  
  reduced_bg = bg_coding[sample(1:length(bg_coding), size = ssize, replace = F, prob = probs)]
  
  return(reduced_bg)
}

function_to_match_CG = function(fg, bg, grain = 0.1, ssize = 2000){
  
  fgseqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, fg)
  bgseqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, bg)
  #we will fit the rpkm of the foreground seqs to a lognormal dist
  #dis$estimate[[1]] is the mean and [[2]] is the sd
  fg_GC_frac = letterFrequency(fgseqs, letters = c("AT", "GC"), as.prob = T)[,2]
  bg_GC_frac = letterFrequency(bgseqs, letters = c("AT", "GC"), as.prob = T)[,2]
  
  dis = fitdist(fg_GC_frac, distr = "lnorm", method = "mme")
  lmean = dis$estimate[[1]]
  lsd = dis$estimate[[2]]
  #pseudocount
  pc = 0.001
  #construct probability vector for sampling
  #use the desnisty function to estimate the probability of each rpkm value
  #each rpkm value in a window of 'grain'
  probs = plnorm(bg_GC_frac+pc+grain, meanlog = lmean, sdlog = lsd) - 
    plnorm(bg_GC_frac+pc, meanlog = lmean, sdlog = lsd)
  
  reduced_bg = bg[sample(1:length(bg), size = ssize, replace = F, prob = probs)]
  
  return(reduced_bg)
}

savePlots = function(filename, fg, bg1, bg_expr, bg_GC){
  
  fgseqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, fg)
  bgseqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, bg1)
  bg_GC_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, bg_GC)
  #we will fit the rpkm of the foreground seqs to a lognormal dist
  #dis$estimate[[1]] is the mean and [[2]] is the sd
  fg_GC_frac = letterFrequency(fgseqs, letters = c("AT", "GC"), as.prob = T)[,2]
  bg_GC_frac = letterFrequency(bgseqs, letters = c("AT", "GC"), as.prob = T)[,2]
  bg_GC_GC_frac = letterFrequency(bg_GC_seqs, letters = c("AT", "GC"), as.prob = T)[,2]
  
  pdf(filename)
  #remember rpkm is V7 in fg and V5 in bg
  pc = 0.01
  ##plot of expression of two groups
  par(mfrow=c(2,1), mar=c(2,2,2,2))
  hist(log(bg1$V5 + pc), breaks = 100, xlab = "rpkm", main = "expression of regions")
  hist(log(fg$V7 + pc), breaks = 100, add = T, col = rgb(1,0,0,1/2))
  legend("topleft", c(paste("bg", ", n=", length(bg1), sep = "") , 
                      paste("fg", ", n=", length(fg), sep = "")), col=c("grey", "red"), lwd=4)  
  hist(log(bg_expr$V5 + pc), breaks = 100, xlab = "rpkm",
       main = "expression of regions, post-filtering", col = rgb(0,0,1,1/2))
  hist(log(fg$V7 + pc), breaks = 100, add = T, col = rgb(1,0,0,1/2))
  legend("topleft", c(paste("bg_expr_control", ", n=", length(bg_expr), sep = "") , 
                      paste("fg", ", n=", length(fg), sep = "")), col=c("blue", "red"), lwd=4)
  
  bg_expr_GC_frac = letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, bg_expr), 
                                    letters = c("AT", "GC"), as.prob = T)[,2]
  
  #plot of GC of two groups
  par(mfrow=c(2,1), mar=c(2,2,2,2))
  hist(bg_expr_GC_frac, breaks = 100, xlim = c(0,1), xlab = "GC content", 
       main = "GC content of regions")
  hist(fg_GC_frac, breaks = 100, add = T, col = rgb(1,0,0,1/2))
  legend("topleft", c(paste("bg_expr_control", ", n=", length(bg_expr), sep = "") , 
                      paste("fg", ", n=", length(fg), sep = "")), col=c("grey", "red"), lwd=4)  
  
  hist(bg_GC_GC_frac, breaks = 100, xlim = c(0,1), xlab = "GC content", 
       main = "GC content of regions, post-filtering", col = rgb(0,0,1,1/2))
  hist(fg_GC_frac, breaks = 100, add = T, col = rgb(1,0,0,1/2))
  legend("topleft", c(paste("bg_expr_match_GC", ", n=", length(bg_GC_GC_frac), sep = "") , 
                      paste("fg", ", n=", length(fg), sep = "")), col=c("blue", "red"), lwd=4)  
  
  #end
  dev.off()
}

ttome_augment = function(RNA_ttome, Pol_ttome, NDR, ext_for_olap = 500){
  
  print(paste0("considering the first this many bp: ", ext_for_olap))
  start.time <- Sys.time()
  RNA5p = promoters(RNA_ttome, upstream = 0, downstream = ext_for_olap)
  Pol5p = promoters(Pol_ttome, upstream = 0, downstream = ext_for_olap)
  pooled = c(RNA5p, Pol5p)
  
  test = mergeByOverlaps(NDR, pooled)
  
  #reset the TUs to their original length
  pooled_unadj = c(RNA_ttome, Pol_ttome)
  names(pooled_unadj) = pooled_unadj$txname
  test$pooled = pooled_unadj[test$pooled$txname]
  
  #one sapply to index into the GRanges for each NDR
  #second sapply to call parse_initiators
  NDR_TUs = sapply(sapply(
    unique(test$id), function(x){test[test$id == x,]$pooled}
  ), parse_initiators)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  return(unlist(GRangesList(NDR_TUs)))
  
}

parse_initiators = function(ranges){
  
  #first, split by strand
  r_plus = ranges[strand(ranges) == "+"]
  r_minus = ranges[strand(ranges) == "-"]
  
  #if one or zero
  if(length(r_plus) < 2){
    r_plus_final = r_plus
  }else{
    ##get self-overlap
    solap = as.table(findOverlaps(r_plus))
    if(sum(solap == 1) > 0){ # at least one singleton
      
      if(all(as.table(findOverlaps(r_plus)) == 1)){    # non-overlapping singetons
        
        sorted = sort(r_plus)
        #give them the same txname... use txname of 'last tx'
        txname_grp = sorted[length(sorted)]$txname
        
        #assign alt to the first one
        alternate = sorted[1]; alternate$alt = T
        r_plus_final = c(alternate, sorted[2:length(sorted)])
        r_plus_final$txname = txname_grp
        
      }else{
        main_TSS = reduce(r_plus[solap > 1])
        #mcols will be lost so need to reassign
        mcols(main_TSS) = as.data.frame(r_plus[solap > 1]) %>% 
          summarise(txname = paste(txname, collapse=';'), 
                    rpkm = mean(rpkm, na.rm=TRUE), 
                    basemeans = mean(basemeans, na.rm=TRUE), 
                    log2FC = mean(log2FC, na.rm=TRUE), 
                    pval = mean(pval, na.rm=TRUE), 
                    Pol_density = mean(Pol_density, na.rm=TRUE))
        
        other_TSS = r_plus[solap == 1]
        #tag it - if it's from the same NDR and same strand, I consider it same tx
        
        sorted = sort(c(other_TSS, main_TSS))
        txname_grp = sorted[length(sorted)]$txname
        ##the first one will be tagged as 'alt'
        alternate = sorted[1]; alternate$alt = T
        r_plus_final = c(alternate, sorted[2:length(sorted)])
        r_plus_final$txname = txname_grp
        
      }
      
    }else{  ##if all overlapping
      
      r_plus_final = reduce(r_plus)
      
      mcols(r_plus_final) = as.data.frame(r_plus) %>% 
        summarise(txname = paste(txname, collapse=';'), 
                  rpkm = mean(rpkm, na.rm=TRUE), 
                  basemeans = mean(basemeans, na.rm=TRUE), 
                  log2FC = mean(log2FC, na.rm=TRUE), 
                  pval = mean(pval, na.rm=TRUE), 
                  Pol_density = mean(Pol_density, na.rm=TRUE))
    }
  }
  
  #now minus strand... this should be a subroutine
  if(length(r_minus) < 2){
    r_minus_final = r_minus
  }else{
    ##get self-overlap
    solap = as.table(findOverlaps(r_minus))
    if(sum(solap == 1) > 0){ # at least one singleton
      
      if(all(as.table(findOverlaps(r_minus)) == 1)){    # non-overlapping singetons
        
        sorted = rev(sort(r_minus))
        #give them the same txname... use txname of 'last tx'
        txname_grp = sorted[length(sorted)]$txname
        
        #assign alt to the first one
        alternate = sorted[1]; alternate$alt = T
        r_minus_final = c(alternate, sorted[2:length(sorted)])
        r_minus_final$txname = txname_grp
        
      }else{
        main_TSS = reduce(r_minus[solap > 1])
        #mcols will be lost so need to reassign
        mcols(main_TSS) = as.data.frame(r_minus[solap > 1]) %>% 
          summarise(txname = paste(txname, collapse=';'), 
                    rpkm = mean(rpkm, na.rm=TRUE), 
                    basemeans = mean(basemeans, na.rm=TRUE), 
                    log2FC = mean(log2FC, na.rm=TRUE), 
                    pval = mean(pval, na.rm=TRUE), 
                    Pol_density = mean(Pol_density, na.rm=TRUE))
        
        other_TSS = r_minus[solap == 1]
        #tag it - if it's from the same NDR and same strand, I consider it same tx
        
        sorted = rev(sort(c(other_TSS, main_TSS)))
        txname_grp = sorted[length(sorted)]$txname
        ##the first one will be tagged as 'alt'
        alternate = sorted[1]; alternate$alt = T
        r_minus_final = c(alternate, sorted[2:length(sorted)])
        r_minus_final$txname = txname_grp
      }
      
    }else{  ##if all overlapping
      
      r_minus_final = reduce(r_minus)
      
      mcols(r_minus_final) = as.data.frame(r_minus) %>% 
        summarise(txname = paste(txname, collapse=';'), 
                  rpkm = mean(rpkm, na.rm=TRUE), 
                  basemeans = mean(basemeans, na.rm=TRUE), 
                  log2FC = mean(log2FC, na.rm=TRUE), 
                  pval = mean(pval, na.rm=TRUE), 
                  Pol_density = mean(Pol_density, na.rm=TRUE))
    }
  }
  
  return(c(r_plus_final, r_minus_final))
  
}

get_5p_GRO = function(ohler, NDR_regions, region, res = 10, filter = 10, exact = T){
  
  #given a region, find TSS based on 5' GROseq data
  # retrieve TSS using NDR
  prom = NDR_regions[NDR_regions$id == region$NDR]
  strand(prom) = strand(region)
  
  candidates = subsetByOverlaps(ohler, prom)
  medval = median(candidates$V5)
  #first filtering
  candidates = candidates[candidates$V5 > medval]
  
  candidates$cluster = csaw::mergeWindows(candidates, ignore.strand = F, tol = res)$ids
  df = as.data.frame(candidates) %>% 
    group_by(cluster) %>% 
    mutate(cluscov = sum(V5)) %>% 
    dplyr::filter(cluscov > filter)
  
  if(nrow(df) == 0){
    #region$ohler = NA
    #make empty ranges
    region$ohler = GRangesList(GRanges())
    return(region)
  }
  #dev.off()
  #hist(unlist(sapply(1:length(candidates), function(x){rep(start(candidates[x]), candidates[x]$V5)})),
  #     breaks = 100)
  
  if(exact == T){
    #get the one with most cov
    final_TSS = makeGRangesFromDataFrame(df %>% 
                                           arrange(desc(V5), .by_group = T) %>% 
                                           dplyr::filter(row_number() == 1) %>% 
                                           ungroup() %>% 
                                           dplyr::select(-c(V4, cluster)), 
                                         keep.extra.columns = T)
    #append to mcols
    region$ohler = GRangesList(final_TSS)
    
    return(region)
  }else{
    ##midpoint
    final_TSS = makeGRangesFromDataFrame(df %>% 
                                           mutate(clus_start = min(start), clus_end = max(end)) %>% 
                                           mutate(clus_width = abs(clus_start - clus_end)) %>% 
                                           arrange(desc(V5), .by_group = T) %>% 
                                           dplyr::filter(row_number() == 1) %>% 
                                           ungroup() %>% 
                                           dplyr::select(-c(V4, V5, cluster, start, end, width)),
                                         start.field = "clus_start",
                                         end.field = "clus_end",
                                         keep.extra.columns = T)
    
    region$ohler = GRangesList(final_TSS)
    
    return(region)
  }
  
}

subroutine_x = function(ranges){
  
  #x = start(ranges)-start(ranges)[1]
  #y = ranges$V5
  data = unlist(sapply(1:length(ranges), function(x){rep(start(ranges[x]), ranges[x]$V5)}))
  par(mfrow = c(1,2))
  hist(data, breaks = 100)
  #repnormmixmodel.sel(do.call('rbind', list(data, data)), k = 1:4)
  #test = normalmixEM(data, lambda = NULL, mu = NULL, sigma = NULL, k = 2, arbvar = F)
  plot(ecdf(data))
  
}

streme_wrapper = function(fg_regions, bg_regions, 
                          minwidth = 4, maxwidth = 6, n_mots = 30, 
                          order = 2, write_dir = "auto",
                          as.bed = T){
  
  #if provided as regions and not sequences
  if(as.bed == T){
    fg_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, fg_regions)
    bg_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, bg_regions)
  }else{
    fg_seqs = fg_regions
    bg_seqs = bg_regions
  }
  
  results = runStreme(
    fg_seqs,
    bg_seqs,
    outdir = write_dir,
    objfun = "de",
    alph = "rna",
    silent = T,
    rna = T,
    minw = minwidth,
    maxw = maxwidth,
    nmotifs = n_mots,
    order = order
  )
  #motlist = results %>% to_list()
  #ggseqlogo::ggseqlogo(lapply(1:n_mots, function(x){motlist[[x]]["motif"]}), ncol = 2)
  
  return(results)
  
}

meme_wrapper = function(fg_regions, bg_regions, 
                        minwidth = 4, maxwidth = 6, n_mots = 30, 
                        order = 2, write_dir = "auto",
                        as.bed = T){
  
  #if provided as regions and not sequences
  if(as.bed == T){
    fg_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, fg_regions)
    bg_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, bg_regions)
  }else{
    fg_seqs = fg_regions
    bg_seqs = bg_regions
  }
  
  results = runMeme(
    fg_seqs,
    neg = bg_seqs,
    outdir = write_dir,
    objfun = "de",
    rna = T,
    minw = minwidth,
    maxw = maxwidth,
    nmotifs = n_mots,
    shuf = order,
    mod = "anr"
  )
  #motlist = results %>% to_list()
  #ggseqlogo::ggseqlogo(lapply(1:n_mots, function(x){motlist[[x]]["motif"]}), ncol = 2)
  
  return(results)
  
}

viewTSS = function(systemTUs_obj, width = 100, shift_val = -50, unique = T, tss_mode = "main",
                   plot_title = NULL){
  
  if(unique == T){
    ohlers = get_uniq_TSS(systemTUs_obj, mode = tss_mode)
    
    consmat = consensusMatrix(getSeq(
      genome, shiftStranded(resize(
        ohlers, 
        width = width, fix = "start"), shift_val)
    ), as.prob = T)[1:4,]
    plot1 = as.data.frame(t(consmat)) %>% 
      mutate(index = row_number()) %>% 
      ggplot(aes(x = index)) + 
      geom_line(aes(y = A), color = "green") + 
      geom_line(aes(y = C), color = "blue") + 
      geom_line(aes(y = G), color = "orange") + 
      geom_line(aes(y = T), color = "red") + 
      xlab("") + ylab("") + ylim(c(0,0.75)) + 
      ggtitle(paste0(plot_title, " (n=", length(ohlers), ")")) + theme_classic()
    print(paste0("running get_main_TSS() to get Ohler TSS w ", tss_mode, " coverage per tx"))
    plot2 = view_motifs(consmat, show.positions = F)
    
    grid.arrange(plot1, plot2)
    
  }else{
    tsses = unlist(GRangesList(systemTUs_obj$ohler))
    
    view_motifs(consensusMatrix(getSeq(
      genome, strandawareshift(resize(tsses, width = width, fix = "start"), shift_val)
    ))[1:4,], show.positions = F)
  }
}

get_uniq_TSS = function(x, mode = "main"){
  
  #given ranges, get the Ohler TSS with highest coverage - ZERO or ONE per tx
  expanded = unlist(GRangesList(x$ohler))
  expanded$tag = names(expanded); names(expanded) = NULL
  
  if(mode == "main"){
    reduced = makeGRangesFromDataFrame(as.data.frame(expanded) %>% 
                                         group_by(tag) %>%
                                         arrange(desc(V5), .by_group = T) %>%
                                         dplyr::filter(row_number() == 1),
                                       keep.extra.columns = T)
  }else if(mode == "leftmost"){
    #for plus strand leftmost, for minus strand rightmost
    exp_plus = expanded[strand(expanded) == "+"]
    exp_minus = expanded[strand(expanded) == "-"]
    exp_plus_reduced = makeGRangesFromDataFrame(as.data.frame(exp_plus) %>% 
                                                  group_by(tag) %>% 
                                                  arrange(start, .by_group = T) %>% 
                                                  dplyr::filter(row_number() == 1),
                                                keep.extra.columns = T)
    exp_minus_reduced = makeGRangesFromDataFrame(as.data.frame(exp_minus) %>% 
                                                   group_by(tag) %>% 
                                                   arrange(desc(end), .by_group = T) %>% 
                                                   dplyr::filter(row_number() == 1),
                                                 keep.extra.columns = T)
    reduced = c(exp_plus_reduced, exp_minus_reduced)
    
  }else if(mode == "rightmost"){
    #for plus strand rightmost, for minus strand leftmost
    exp_plus = expanded[strand(expanded) == "+"]
    exp_minus = expanded[strand(expanded) == "-"]
    exp_plus_reduced = makeGRangesFromDataFrame(as.data.frame(exp_plus) %>% 
                                                  group_by(tag) %>% 
                                                  arrange(desc(end), .by_group = T) %>% 
                                                  dplyr::filter(row_number() == 1),
                                                keep.extra.columns = T)
    exp_minus_reduced = makeGRangesFromDataFrame(as.data.frame(exp_minus) %>% 
                                                   group_by(tag) %>% 
                                                   arrange(start, .by_group = T) %>% 
                                                   dplyr::filter(row_number() == 1),
                                                 keep.extra.columns = T)
    reduced = c(exp_plus_reduced, exp_minus_reduced)
    
  }else if(mode == "lowest"){
    reduced = makeGRangesFromDataFrame(as.data.frame(expanded) %>% 
                                         group_by(tag) %>%
                                         arrange(desc(V5), .by_group = T) %>%
                                         dplyr::filter(row_number() == n()),
                                       keep.extra.columns = T)
  }
  print(paste0(length(expanded), " initial TSS, ", length(reduced), " final TSS"))
  
  return(reduced)
}

CLIP_systemTUs_get = function(clip_peaks, systemtus){
  
  #return systemTUs object with associated CLIP peaks and SCORE as metadata
  clop = findOverlaps(systemtus, clip_peaks, maxgap = 1000, ignore.strand = F)
  t1 = systemtus[queryHits(clop)]
  t2 = clip_peaks[subjectHits(clop)]
  t2$NDR_tag = t1$NDR_tag
  
  summarised = as.data.frame(t2) %>% 
    group_by(NDR_tag) %>% 
    summarise(mean_adj_score = mean(adj_score), mean_logFC = mean(log2FC))
  
  CLIP_TUs = systemtus[summarised$NDR_tag]
  CLIP_TUs$mean_adj_score = summarised$mean_adj_score
  CLIP_TUs$mean_log2FC = summarised$mean_logFC
  
  
  return(CLIP_TUs[order(CLIP_TUs$mean_adj_score, decreasing = T)])
}

get_distance_clip = function(systemtus, clipper_peaks){
  ###now, for the clipper peaks, record distance to nearest TSS
  ##initialise metadata columns
  clipper_peaks$proseq_tss = GRangesList(GRanges())
  clipper_peaks$ohler_tss = GRangesList(GRanges())
  clipper_peaks$NDR_tag = NA
  
  #find tx it maps to
  olap = findOverlaps(systemtus, clipper_peaks, maxgap = 1000, ignore.strand = F)
  
  clipper_peaks[subjectHits(olap)[!duplicated(subjectHits(olap))]]$NDR_tag = 
    systemtus[queryHits(olap)[!duplicated(subjectHits(olap))]]$NDR_tag
  
  clipper_peaks[!is.na(clipper_peaks$NDR_tag)]$ohler_tss = 
    systemtus[clipper_peaks$NDR_tag[!is.na(clipper_peaks$NDR_tag)]]$ohler
  
  proseq_starts = resize(systemtus[clipper_peaks$NDR_tag[!is.na(clipper_peaks$NDR_tag)]], 
                         width = 1, fix = "start")
  mcols(proseq_starts) = NULL
  
  clipper_peaks[!is.na(clipper_peaks$NDR_tag)]$proseq_tss = proseq_starts
  
  ###distance calculation
  dist_to_ohler = lapply(1:length(clipper_peaks), function(x){
    distanceToNearest(clipper_peaks[x], clipper_peaks[x]$ohler_tss[[1]], ignore.strand=FALSE)
  })
  
  dist_to_proseq = lapply(1:length(clipper_peaks), function(x){
    distanceToNearest(clipper_peaks[x], clipper_peaks[x]$proseq_tss[[1]], ignore.strand=FALSE)
  })
  
  clipper_peaks$dist2ohler = dist_to_ohler
  clipper_peaks$dist2pro = dist_to_proseq
  ###return peaks with the additional info
  
  return(clipper_peaks)
  
}

## function to output TSS seqlogo and do motif search using streme
TSS_and_motif = function(outdir, filename, fg, bg, shift_val, width, fg_tss_mode = "main", bg_tss_mode = "main",
                         do_heatmap = F, usePROSEQ = F, minw = 4, maxw = 10, streme_order = 1,
                         fg_name = NULL, bg_name = NULL, usestreme = T){
  
  ##make the outdir
  dir.create(file.path(getwd(), "Shannahan", outdir))
  
  pdf(file.path(getwd(), "Shannahan", outdir, filename), height = 5, width = 10)
  ###viewTSS
  p1 = viewTSS(fg, width = width, shift_val = shift_val, unique = T, 
               tss_mode = fg_tss_mode, plot_title = fg_name)
  print(p1)
  p2 = viewTSS(bg, width = width, shift_val = shift_val, unique = T, 
               tss_mode = bg_tss_mode, plot_title = bg_name)
  print(p2)
  
  ##heatmap
  if(do_heatmap == T){
    
    print("making heatmap of regions")
    testseqs = getSeq(
      genome, shiftStranded(resize(
        get_uniq_TSS(fg, mode = tss_mode), 
        width = width, fix = "start"), shift_val)
    )
    testmat = as.matrix(testseqs)
    
    testmat[testmat == "-"] = 0
    testmat[testmat == "A"] = 1
    testmat[testmat == "C"] = 2
    testmat[testmat == "G"] = 2
    testmat[testmat == "T"] = 1
    
    testmat = apply(testmat, 2, as.numeric)
    
    pheatmap(testmat, cluster_rows = F, cluster_cols = F,
             color = c("white", "green", "blue", "orange", "red"),
             breaks = c(0, 0.5, 1.5, 2.5, 3.5, 4.5))  # distances 0 to 0.5 are white, etc...
    
    testseqs2 = getSeq(
      genome, shiftStranded(resize(
        get_uniq_TSS(bg, mode = tss_mode), 
        width = width, fix = "start"), shift_val)
    )
    testmat2 = as.matrix(testseqs2)
    
    testmat2[testmat2 == "-"] = 0
    testmat2[testmat2 == "A"] = 1
    testmat2[testmat2 == "C"] = 2
    testmat2[testmat2 == "G"] = 2
    testmat2[testmat2 == "T"] = 1
    
    testmat2 = apply(testmat2, 2, as.numeric)
    
    pheatmap(testmat2, cluster_rows = F, cluster_cols = F,
             color = c("white", "green", "blue", "orange", "red"),
             breaks = c(0, 0.5, 1.5, 2.5, 3.5, 4.5))
    dev.off()
  }else{
    dev.off()
  }
  
  if(usePROSEQ == T){
    print("using proseq start sites for motif search")
    fg_seqs = getSeq(genome, promoters(fg, upstream = 0, downstream = width))
    bg_seqs = getSeq(genome, promoters(bg, upstream = 0, downstream = width))
  }else{
    fg_seqs = getSeq(
      genome, shiftStranded(resize(
        get_uniq_TSS(fg, mode = fg_tss_mode), 
        width = width, fix = "start"), shift_val)
    )
    names(fg_seqs) = get_uniq_TSS(fg, mode = fg_tss_mode)$tag
    bg_seqs = getSeq(
      genome, shiftStranded(resize(
        get_uniq_TSS(bg, mode = bg_tss_mode), 
        width = width, fix = "start"), shift_val)
    )
    names(bg_seqs) = get_uniq_TSS(bg, mode = bg_tss_mode)$tag
  }
  
  ##motif search
  if(usestreme == T){
    print("running streme")
    
    sapply(streme_order, function(x){
      streme_wrapper(fg_seqs, 
                     bg_seqs, 
                     write_dir = paste0(file.path(getwd(), "Shannahan", outdir, "streme_order_"),
                                        as.character(x)),
                     order = x,
                     minwidth = minw,
                     maxwidth = maxw,
                     as.bed = F)
    })
  }else{
    print("running meme")
    
    sapply(streme_order, function(x){
      meme_wrapper(fg_seqs, 
                   bg_seqs, 
                   write_dir = paste0(file.path(getwd(), "Shannahan", outdir, "meme_order_"),
                                      as.character(x)),
                   order = x,
                   minwidth = minw,
                   maxwidth = maxw,
                   as.bed = F)
    })
  }
  
  
  
}

get_proseq_tss = function(systemtus, proseq_macspeaks){
  
  ###NDR based
  #for now, take the highest wscoring proseq tss
  systemtus_w_proseq = lapply(1:length(systemtus), function(x){
    parse_one_proseq_tss(proseq_macspeaks, systemtus[x])
  })
  
  return(systemtus_w_proseq)
}

parse_one_proseq_tss = function(proseq_macspeaks, one_tu){
  
  #first check if exists
  if(one_tu$NDR %in% proseq_macspeaks$NDR == F){
    one_tu$proseq_tss = GRangesList(GRanges())
    return(one_tu)
  }
  proseq_region = proseq_macspeaks[proseq_macspeaks$NDR == one_tu$NDR & strand(proseq_macspeaks) == strand(one_tu)]
  if(length(proseq_region) == 0){
    one_tu$proseq_tss = GRangesList(GRanges())
    return(one_tu)
  }
  #in case there are more
  tss_reg = proseq_region[order(proseq_region$pval, decreasing = T)][1]
  #tss = promoters(tss_reg, upstream = 0, downstream = 1)
  one_tu$proseq_tss = GRangesList(tss_reg)
  return(one_tu)
}

writeregions = function(systemtus_list, ups, downs, tags,
                        dir = "/Users/ieo5559/enhancedClip/Shannahan/region_fastas/",
                        PROSEQ = F){
  
  stopifnot(length(systemtus_list) == length(tags))
  
  if(PROSEQ == F){
    reglist = lapply(1:length(systemtus_list), function(x){
      promoters(get_uniq_TSS(systemtus_list[[x]], mode = "main"), 
                upstream = ups, downstream = downs)
    })
    
    for(i in 1:length(reglist)){
      names(reglist[[i]]) = reglist[[i]]$tag
    }
  }else{
    reglist = lapply(1:length(systemtus_list), function(x){
      promoters(systemtus_list[[x]], upstream = ups, downs = downs)
    })
    
    for(i in 1:length(reglist)){
      names(reglist[[i]]) = reglist[[i]]$NDR_tag
    }
  }
  
  
  for(j in 1:length(reglist)){
    writeXStringSet(
      getSeq(genome, reglist[[j]]), paste0(dir, tags[[j]], ".fa")
    )
    write.table(as.data.frame(reglist[[j]])[,c(1,2,3,4,6,5)], 
                file = paste0(dir, tags[[j]], ".bed"),
                quote = F, sep = "\t", 
                row.names = F,
                col.names = F)
  }
  
  
}

shiftStranded = function(x, value=0L,...) GenomicRanges::shift(x ,shift=value*ifelse('-'==strand(x),-1,1),...)

RNAcount = function(regions, factor=1e9){
  
  bams_f = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_fwd_clean.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/94_SE_fwd_clean.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/95_SE_fwd_clean.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/07/fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/15/fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/23/fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/08/fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/16/fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/24/fwd.bam")
  
  bams_r = c("/hpcnfs/data/GN2/gmandana/take2project/clip/analysis/read2/SE/70_SE_rev_clean.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/94_SE_rev_clean.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/SE/95_SE_rev_clean.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/07/rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/15/rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/23/rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/08/rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/16/rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/24/rev.bam")
  
  param = readParam(minq=10, pe="both")
  
  if(is.null(names(regions))){
    names(regions) = paste("id", 1:length(regions), sep = "")
  }
  
  regions_plus = regions[strand(regions)=="+"]
  regions_minus = regions[strand(regions)=="-"]
  
  #count for and rev separately - count ONLY RNA SEQ BAMS, check param
  regcounts_plus <- regionCounts(bams_f[4:9], regions_plus, param=param, BPPARAM = MulticoreParam())
  regcounts_minus <- regionCounts(bams_r[4:9], regions_minus, param=param, BPPARAM = MulticoreParam())
  
  #allcounts = rbind(assay(regcounts_plus), assay(regcounts_minus))
  
  cp = sweep(sweep(assay(regcounts_plus), MARGIN = 2, STATS = regcounts_plus$totals, FUN = "/"), MARGIN = 2, STATS = width(regions_plus), FUN = "/")
  cm = sweep(sweep(assay(regcounts_minus), MARGIN = 2, STATS = regcounts_minus$totals, FUN = "/"), MARGIN = 2, STATS = width(regions_minus), FUN = "/")
  
  means_p = cbind(rowMeans(cp[,1:3]), rowMeans(cp[,4:6]))
  names(means_p) = names(regions_plus)
  means_m = cbind(rowMeans(cm[,1:3]), rowMeans(cm[,4:6]))
  names(means_m) = names(regions_minus)
  
  allcounts = rbind(means_p*factor, means_m*factor)
  
  countstable = allcounts[names(regions),]
  
  return(allcounts[names(regions),])
  
}

fig_2f_generation = function(clip_tss, bg_tss, outfilename, rbns_pwm, fimobgfile){
  
  regions_fg = promoters(clip_tss, upstream = 0, downstream = 100)
  regions_bg = promoters(bg_tss, upstream = 0, downstream = 100)
  
  fimo_bg_generator(regions_fg, fimobgfile)
  thresh_list = seq(1e-3, 5e-2, length.out = 100)
  pdf(outfilename)
  
  output_m1 = lapply(thresh_list, function(x){
    length(unique(as.data.frame(runFimo(sequences = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, regions_fg),
                                   motifs = rbns_pwm[1], norc = T, thresh = x, bfile = fimobgfile,
                                   max_stored_scores = 100000,
                                   parse_genomic_coord = F,
                                   silent = T))$seqnames))
  })

  output_m2 = lapply(thresh_list, function(x){
    length(unique(as.data.frame(runFimo(sequences = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, regions_fg),
                                        motifs = rbns_pwm[2], norc = T, thresh = x, bfile = fimobgfile,
                                        max_stored_scores = 100000,
                                        parse_genomic_coord = F,
                                        silent = T))$seqnames))
  })

  output_both = lapply(thresh_list, function(x){
    length(unique(as.data.frame(runFimo(sequences = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, regions_fg),
                                        motifs = c(rbns_pwm[1], rbns_pwm[2]), norc = T, thresh = x, bfile = fimobgfile,
                                        max_stored_scores = 100000,
                                        parse_genomic_coord = F,
                                        silent = T))$seqnames))
  })


  fimo_matrix = matrix(c(unlist(output_m1), unlist(output_m2), unlist(output_both)), ncol = 3, byrow = F)

  #matplot(thresh_list, fimo_matrix, type = "l", ylim = c(0,length(clip_proseq_tss)))

  ##ggplot2

  mot_plot = as.data.frame(fimo_matrix) %>%
    `colnames<-`(c("motif1", "motif2", "motif1and2")) %>%
    melt() %>%
    mutate(thresh_vals = rep(thresh_list, 3)) %>%
    ggplot(aes(x=thresh_vals, y=value, group=variable, color=variable)) +
    geom_line() +
    ylab("number of unique peak sequences with motif") +
    xlab("FIMO p-value threshold") +
    theme_classic()

  print(mot_plot)

  ## now the other one
  
  output_100 = lapply(thresh_list, function(x){
    length(unique(seqnames(runFimo(sequences = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, 
                                                           promoters(clip_tss, upstream = 0, downstream = 100)), 
                                        motifs = c(rbns_pwm[1], rbns_pwm[2]), norc = T, thresh = x, bfile = fimobgfile,
                                        max_stored_scores = 100000,
                                        parse_genomic_coord = F,
                                        silent = T))))
  })
  
  fimo_bg_generator(promoters(clip_tss, upstream = 0, downstream = 150), fimobgfile)
  
  output_150 = lapply(thresh_list, function(x){
    length(unique(as.data.frame(runFimo(sequences = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, 
                                                           promoters(clip_tss, upstream = 0, downstream = 150)), 
                                        motifs = c(rbns_pwm[1], rbns_pwm[2]), norc = T, thresh = x, bfile = fimobgfile,
                                        max_stored_scores = 100000,
                                        parse_genomic_coord = F,
                                        silent = T))$seqnames))
  })
  
  fimo_bg_generator(promoters(clip_tss, upstream = 0, downstream = 200), fimobgfile)
  
  output_200 =  lapply(thresh_list, function(x){
    length(unique(as.data.frame(runFimo(sequences = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, 
                                                           promoters(clip_tss, upstream = 0, downstream = 200)), 
                                        motifs = c(rbns_pwm[1], rbns_pwm[2]), norc = T, thresh = x, bfile = fimobgfile,
                                        max_stored_scores = 100000,
                                        parse_genomic_coord = F,
                                        silent = T))$seqnames))
  })
  
  extension_matrix = matrix(c(unlist(output_100), unlist(output_150), unlist(output_200)), ncol = 3, byrow = F)
  
  ext_plot = as.data.frame(extension_matrix) %>% 
    `colnames<-`(c("100bp", "150bp", "200bp")) %>% 
    melt() %>%
    mutate(thresh_vals = rep(thresh_list, 3)) %>% 
    ggplot(aes(x=thresh_vals, y=value, group=variable, color=variable)) + 
    geom_line() + 
    ylab("number of unique peak sequences with motif") + 
    xlab("FIMO p-value threshold") +
    theme_classic()
  
  print(ext_plot)
  
  dev.off()
  
}

# pass list for bwfiles (plus and minus as vectors)
getcovmat = function(bwfiles, regions, nbins, norm = T, smoothwin = 10){
  
  if(is.null(names(regions))){
    names(regions) = paste0("id", 1:length(regions))
  }
  
  if(length(unique(width(regions))) != 1){
    stop("need to be same length")
  }
  
  samnames = basename(file_path_sans_ext(bwfiles[[1]]))
  
  r_p = regions[strand(regions) == "+"]
  r_m = regions[strand(regions) == "-"]
  
  plusmat = do.call('rbind', lapply(1:length(bwfiles[[1]]), function(x){
    colMeans(bigwig_collect(bwfiles[[1]][x], regions = r_p, norm = norm, nbins = nbins),
             na.rm = T)
  }))
  
  minmat = do.call('rbind', lapply(1:length(bwfiles[[2]]), function(x){
    colMeans(bigwig_collect(bwfiles[[2]][x], regions = r_m, norm = norm, nbins = nbins),
             na.rm = T)
  }))
  
  covmat = (minmat + plusmat)/2
  #final_matrix = t(covmat)
  rownames(covmat) = samnames
  
  return(covmat)
  
}

#wrapper for a wrapper....
regionscounts = function(bams, regions, parameter = readParam(minq=10, pe="both")){
  
  #split regions
  if(is.null(regions$id)){
    regions$id = paste("id", 1:length(regions), sep = "")
  }
  
  regions_plus = regions[strand(regions)=="+"]
  regions_minus = regions[strand(regions)=="-"]
  
  regcounts_plus <- regionCounts(bams[[1]], regions_plus, param=parameter, BPPARAM = MulticoreParam())
  regcounts_minus <- regionCounts(bams[[2]], regions_minus, param=parameter, BPPARAM = MulticoreParam())
  
  rownames(regcounts_plus) = regions_plus$id
  rownames(regcounts_minus) = regions_minus$id
  
  tot = rbind(assay(regcounts_plus), assay(regcounts_minus))
  
  return(list("counts" = tot[regions$id,], "totals" = regcounts_plus$totals + regcounts_minus$totals))
  
}

shiftStranded = function(x, value=0L,...) GenomicRanges::shift(x ,shift=value*ifelse('-'==strand(x),-1,1),...)


CLIP_peak_clustering_v1 = function(peaks, 
                                   siWfile = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/T3_ghostbusters/siW__better_annotated_20Sep.bed"){
  
  bams_f = c("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/07/fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/15/fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/23/fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/08/fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/16/fwd.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/24/fwd.bam")
  
  bams_r = c("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/07/rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/15/rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/23/rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/08/rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/16/rev.bam",
             "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/mapped/24/rev.bam")
  
  siW = customimportbed6(siWfile)
  #colnames(mcols(siW)) = c("pval", "FDR", "logFC", "ctl_rpkm")
  names(mcols(siW)) = c("name", "ctl_rpkm", "pval", "FDR", "logFC")
  siW$id = paste0("id", 1:length(siW))
  
  clusters = as.data.frame(peaks) %>% 
    group_by(txname) %>% 
    mutate(cluster_start = min(start)) %>% 
    mutate(cluster_end = max(end)) %>% 
    mutate(coords = paste(seqnames, paste(start, end, sep = "-"), sep = ":")) %>% 
    dplyr::select(-c("start", "end", "width")) %>% 
    summarise_all(paste, collapse = ",")
  
  #for later
  mean_scores = as.data.frame(peaks) %>% 
    group_by(txname) %>% 
    summarise(cluster_mean_score = mean(adj_score))
  
  clusters$seqnames = unlist(lapply(clusters$seqnames, 
                                    function(x){str_split(x, ",")[[1]][[1]]}))
  clusters$strand = unlist(lapply(clusters$strand, 
                                  function(x){str_split(x, ",")[[1]][[1]]}))
  clusters$cluster_start = as.numeric(unlist(lapply(clusters$cluster_start, 
                                                    function(x){str_split(x, ",")[[1]][[1]]})))
  clusters$cluster_end = as.numeric(unlist(lapply(clusters$cluster_end, 
                                                  function(x){str_split(x, ",")[[1]][[1]]})))
  
  clusters = clusters[,c(2,8,9,1,3,4,5,6,7,10)]
  
  clusters_ranges = makeGRangesFromDataFrame(clusters,
                                             start.field = "cluster_start",
                                             end.field = "cluster_end",
                                             strand.field = "strand",
                                             keep.extra.columns = T)
  
  clus_olap = findOverlaps(clusters_ranges, siW, maxgap = 100, ignore.strand = F)
  
  ##append the overlap info
  clusters[queryHits(clus_olap)[!duplicated(queryHits(clus_olap))],]
  
  toappend = as.data.frame(mcols(siW[subjectHits(clus_olap)[!duplicated(queryHits(clus_olap))]])) %>% 
    mutate(summ = paste(name, FDR, id, sep = ";"))
  
  clusters$overlap = NA
  clusters[queryHits(clus_olap)[!duplicated(queryHits(clus_olap))],]$overlap = toappend$summ
  
  clusters = as_tibble(merge(clusters, mean_scores))
  
  ###count at the cluster level
  
  regions_plus = promoters(clusters_ranges[strand(clusters_ranges)=="+"], upstream = 0,
                           downstream = 1000)
  regions_minus = promoters(clusters_ranges[strand(clusters_ranges)=="-"], upstream = 0,
                            downstream = 1000)
  
  #count for and rev separately - count ONLY RNA SEQ BAMS, check param
  param = readParam(minq=10, pe="both")
  regcounts_plus <- regionCounts(bams_f, regions_plus, param=param, BPPARAM = MulticoreParam())
  regcounts_minus <- regionCounts(bams_r, regions_minus, param=param, BPPARAM = MulticoreParam())
  
  plusmat = sweep(x = assay(regcounts_plus), MARGIN = 2, STATS = regcounts_plus$totals, FUN = "/")
  minmat = sweep(x = assay(regcounts_minus), MARGIN = 2, STATS = regcounts_minus$totals, FUN = "/")
  
  regions_plus$ctlcounts = rowMeans(plusmat[,c(1:3)])
  regions_plus$auxcounts = rowMeans(plusmat[,c(4:6)])
  
  regions_minus$ctlcounts = rowMeans(minmat[,c(1:3)])
  regions_minus$auxcounts = rowMeans(minmat[,c(4:6)])
  
  regions_all = c(regions_plus, regions_minus)
  
  clusters = as_tibble(merge(clusters, as.data.frame(regions_all)[,c(6,12,13)]))
  
  clusters_final = clusters[order(clusters$cluster_mean_score, decreasing = T),] %>% 
    mutate(t1 = log2(ctlcounts*1e6), t2 = log2(auxcounts*1e6))
  
  return(clusters_final)
  
  
}



# this function normalises EACH REGION i.e., intensity for ALL SAMPLES 
# for each region are scaled to 1
getcovmat2 = function(bwfiles, regions, nbins, norm = F, smoothwin = 10){
  
  if(is.null(names(regions))){
    names(regions) = paste0("id", 1:length(regions))
  }
  
  if(length(unique(width(regions))) != 1){
    stop("need to be same length")
  }
  
  samnames = basename(file_path_sans_ext(bwfiles[[1]]))
  
  r_p = regions[strand(regions) == "+"]
  r_m = regions[strand(regions) == "-"]
  
  plusmat = lapply(1:length(bwfiles[[1]]), function(x){
    bigwig_collect(bwfiles[[1]][x], regions = r_p, norm = F, nbins = nbins)
  })
  
  plusarray = array(unlist(plusmat), c(length(r_p),nbins,length(bwfiles[[1]])))
  
  #1st region, all bins, 7th sample
  #matplot(testarray[1,,7], type = "l")
  #these are your sizefactors
  SF = sapply(1:length(r_p), function(x){max(plusarray[x,,], na.rm = T)})
  ##sweep
  plusarray_scaled = sweep(plusarray, MARGIN = 1, FUN = "/", STATS = SF)
  ##now collapse the matrix
  plusarray_scaled_coll = apply(plusarray_scaled, c(2,3), FUN = mean, na.rm = T)
  
  ######now same for reverse
  minmat = lapply(1:length(bwfiles[[2]]), function(x){
    bigwig_collect(bwfiles[[2]][x], regions = r_m, norm = F, nbins = nbins)
  })
  
  minarray = array(unlist(minmat), c(length(r_m),nbins,length(bwfiles[[2]])))
  SF_m = sapply(1:length(r_m), function(x){max(minarray[x,,], na.rm = T)})
  minarray_scaled = sweep(minarray, MARGIN = 1, FUN = "/", STATS = SF_m)
  minarray_scaled_coll = apply(minarray_scaled, c(2,3), FUN = mean, na.rm = T)
  
  covmat = (minarray_scaled_coll + plusarray_scaled_coll)/2
  
  final_matrix = t(covmat)
  rownames(final_matrix) = samnames
  
  return(final_matrix)
  
}


##### for regions of differing length
covmat_scalereg = function(bwfiles, regions, numbins){
  
  if(is.null(names(regions))){
    names(regions) = paste0("id", 1:length(regions))
  }
  
  r_p = regions[strand(regions) == "+"]
  r_m = regions[strand(regions) == "-"]
  
  samnames = basename(file_path_sans_ext(bwfiles[[1]]))
  
  plusmat = lapply(1:length(bwfiles[[1]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[1]][x], gr = r_p, Nbin = NULL), 
                   summarize, n = numbins)
    )
  })
  
  plusarray = array(unlist(plusmat), c(length(r_p),numbins,length(bwfiles[[1]])))
  
  #1st region, all bins, 7th sample
  #matplot(testarray[1,,7], type = "l")
  #these are your sizefactors
  SF = sapply(1:length(r_p), function(x){max(plusarray[x,,], na.rm = T)})
  ##sweep
  plusarray_scaled = sweep(plusarray, MARGIN = 1, FUN = "/", STATS = SF)
  ##now collapse the matrix
  plusarray_scaled_coll = apply(plusarray_scaled, c(2,3), FUN = mean, na.rm = T)
  
  ######now same for reverse
  minmat = lapply(1:length(bwfiles[[2]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[2]][x], gr = r_m, Nbin = NULL), 
                   summarize, n = numbins)
    )
  })
  
  minarray = array(unlist(minmat), c(length(r_m), numbins, length(bwfiles[[2]])))
  SF_m = sapply(1:length(r_m), function(x){max(minarray[x,,], na.rm = T)})
  minarray_scaled = sweep(minarray, MARGIN = 1, FUN = "/", STATS = SF_m)
  minarray_scaled_coll = apply(minarray_scaled, c(2,3), FUN = mean, na.rm = T)
  
  minarray_scaled_coll_corr = apply(minarray_scaled_coll, 2, FUN = rev)
  
  covmat = (minarray_scaled_coll_corr + plusarray_scaled_coll)/2
  
  return(t(covmat))
  
}


covmat_scalereg_w_flanks = function(bwfiles, regions, numbins, ups, downs,
                                    upbins, downbins){
  
  if(is.null(names(regions))){
    names(regions) = paste0("id", 1:length(regions))
  }
  
  r_p = regions[strand(regions) == "+"]
  r_m = regions[strand(regions) == "-"]
  
  samnames = basename(file_path_sans_ext(bwfiles[[1]]))
  
  #browser()
  
  plusmat = lapply(1:length(bwfiles[[1]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[1]][x], gr = r_p, Nbin = NULL), 
                   summarize, n = numbins)
    )
  })
  
  plusarray = array(unlist(plusmat), c(length(r_p),numbins,length(bwfiles[[1]])))
  
  #1st region, all bins, 7th sample
  #matplot(testarray[1,,7], type = "l")
  #these are your sizefactors
  SF = sapply(1:length(r_p), function(x){max(plusarray[x,,], na.rm = T)})
  ##sweep
  plusarray_scaled = sweep(plusarray, MARGIN = 1, FUN = "/", STATS = SF)
  ##now collapse the matrix
  plusarray_scaled_coll = apply(plusarray_scaled, c(2,3), FUN = mean, na.rm = T)
  
  ######now same for reverse
  minmat = lapply(1:length(bwfiles[[2]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[2]][x], gr = r_m, Nbin = NULL), 
                   summarize, n = numbins)
    )
  })
  
  minarray = array(unlist(minmat), c(length(r_m), numbins, length(bwfiles[[2]])))
  SF_m = sapply(1:length(r_m), function(x){max(minarray[x,,], na.rm = T)})
  minarray_scaled = sweep(minarray, MARGIN = 1, FUN = "/", STATS = SF_m)
  minarray_scaled_coll = apply(minarray_scaled, c(2,3), FUN = mean, na.rm = T)
  
  minarray_scaled_coll_corr = apply(minarray_scaled_coll, 2, FUN = rev)
  
  covmat = (minarray_scaled_coll_corr + plusarray_scaled_coll)/2
  
  ####flanks
  
  r_p_ups = promoters(r_p, upstream = ups, downstream = 0)
  r_m_ups = promoters(r_m, upstream = ups, downstream = 0)
  
  plusmat_ups = lapply(1:length(bwfiles[[1]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[1]][x], gr = r_p_ups, Nbin = NULL), 
                   summarize, n = upbins)
    )
  })
  
  plusarray_ups = array(unlist(plusmat_ups), c(length(r_p_ups),upbins,length(bwfiles[[1]])))
  plusarray_ups_scaled = sweep(plusarray_ups, MARGIN = 1, FUN = "/", STATS = SF) #prev SF
  plusarray_ups_scaled_coll = apply(plusarray_ups_scaled, c(2,3), FUN = mean, na.rm = T)
  
  minmat_ups = lapply(1:length(bwfiles[[2]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[2]][x], gr = r_m_ups, Nbin = NULL), 
                   summarize, n = upbins)
    )
  })
  
  minarray_ups = array(unlist(minmat_ups), c(length(r_m_ups),upbins,length(bwfiles[[2]])))
  minarray_ups_scaled = sweep(minarray_ups, MARGIN = 1, FUN = "/", STATS = SF_m) #prev SF
  minarray_ups_scaled_coll = apply(apply(minarray_ups_scaled, c(2,3), FUN = mean, na.rm = T),
                                   2, FUN = rev)
  
  ######## downs
  r_p_downs = promoters(resize(r_p, fix = "end", width = 1), upstream = 0, downstream = downs)
  r_m_downs = promoters(resize(r_m, fix = "end", width = 1), upstream = 0, downstream = downs)

  plusmat_downs = lapply(1:length(bwfiles[[1]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[1]][x], gr = r_p_downs, Nbin = NULL), 
                   summarize, n = downbins)
    )
  })
  
  plusarray_downs = array(unlist(plusmat_downs), c(length(r_p_downs),downbins,length(bwfiles[[1]])))
  plusarray_downs_scaled = sweep(plusarray_downs, MARGIN = 1, FUN = "/", STATS = SF) #prev SF
  plusarray_downs_scaled_coll = apply(plusarray_downs_scaled, c(2,3), FUN = mean, na.rm = T)
  
  minmat_downs = lapply(1:length(bwfiles[[2]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[2]][x], gr = r_m_downs, Nbin = NULL), 
                   summarize, n = downbins)
    )
  })
  
  minarray_downs = array(unlist(minmat_downs), c(length(r_m_downs),downbins,length(bwfiles[[2]])))
  minarray_downs_scaled = sweep(minarray_downs, MARGIN = 1, FUN = "/", STATS = SF_m) #prev SF
  minarray_downs_scaled_coll = apply(apply(minarray_downs_scaled, c(2,3), FUN = mean, na.rm = T),
                                   2, FUN = rev)
  
  
  #browser()
  #####stitch
  
  covmat_ups = (minarray_ups_scaled_coll + plusarray_ups_scaled_coll)/2
  covmat = (minarray_scaled_coll_corr + plusarray_scaled_coll)/2
  covmat_downs = (minarray_downs_scaled_coll + plusarray_downs_scaled_coll)/2
  
  final_matrix = t(rbind(covmat_ups, covmat, covmat_downs))
  rownames(final_matrix) = samnames
  
  return(final_matrix)
  
}

covmat_noscalereg_w_flanks = function(bwfiles, regions, numbins, ups, downs,
                                    upbins, downbins){
  
  if(is.null(names(regions))){
    names(regions) = paste0("id", 1:length(regions))
  }
  
  r_p = regions[strand(regions) == "+"]
  r_m = regions[strand(regions) == "-"]
  
  samnames = basename(file_path_sans_ext(bwfiles[[1]]))
  
  #browser()
  
  plusmat = lapply(1:length(bwfiles[[1]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[1]][x], gr = r_p, Nbin = NULL), 
                   summarize, n = numbins)
    )
  })
  
  plusarray = array(unlist(plusmat), c(length(r_p),numbins,length(bwfiles[[1]])))
  
  ##now collapse the matrix
  plusarray_scaled_coll = apply(plusarray, c(2,3), FUN = mean, na.rm = T)
  
  ######now same for reverse
  minmat = lapply(1:length(bwfiles[[2]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[2]][x], gr = r_m, Nbin = NULL), 
                   summarize, n = numbins)
    )
  })
  
  minarray = array(unlist(minmat), c(length(r_m), numbins, length(bwfiles[[2]])))
  minarray_scaled_coll = apply(minarray, c(2,3), FUN = mean, na.rm = T)
  
  minarray_scaled_coll_corr = apply(minarray_scaled_coll, 2, FUN = rev)
  
  covmat = (minarray_scaled_coll_corr + plusarray_scaled_coll)/2
  
  ####flanks
  
  r_p_ups = promoters(r_p, upstream = ups, downstream = 0)
  r_m_ups = promoters(r_m, upstream = ups, downstream = 0)
  
  plusmat_ups = lapply(1:length(bwfiles[[1]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[1]][x], gr = r_p_ups, Nbin = NULL), 
                   summarize, n = upbins)
    )
  })
  
  plusarray_ups = array(unlist(plusmat_ups), c(length(r_p_ups),upbins,length(bwfiles[[1]])))
  
  plusarray_ups_scaled_coll = apply(plusarray_ups, c(2,3), FUN = mean, na.rm = T)
  
  minmat_ups = lapply(1:length(bwfiles[[2]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[2]][x], gr = r_m_ups, Nbin = NULL), 
                   summarize, n = upbins)
    )
  })
  
  minarray_ups = array(unlist(minmat_ups), c(length(r_m_ups),upbins,length(bwfiles[[2]])))
  
  minarray_ups_scaled_coll = apply(apply(minarray_ups, c(2,3), FUN = mean, na.rm = T),
                                   2, FUN = rev)
  
  ######## downs
  r_p_downs = promoters(resize(r_p, fix = "end", width = 1), upstream = 0, downstream = downs)
  r_m_downs = promoters(resize(r_m, fix = "end", width = 1), upstream = 0, downstream = downs)
  
  plusmat_downs = lapply(1:length(bwfiles[[1]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[1]][x], gr = r_p_downs, Nbin = NULL), 
                   summarize, n = downbins)
    )
  })
  
  plusarray_downs = array(unlist(plusmat_downs), c(length(r_p_downs),downbins,length(bwfiles[[1]])))
  
  plusarray_downs_scaled_coll = apply(plusarray_downs, c(2,3), FUN = mean, na.rm = T)
  
  minmat_downs = lapply(1:length(bwfiles[[2]]), function(x){
    do.call('rbind', 
            lapply(coverageFromBW(bwfile = bwfiles[[2]][x], gr = r_m_downs, Nbin = NULL), 
                   summarize, n = downbins)
    )
  })
  
  minarray_downs = array(unlist(minmat_downs), c(length(r_m_downs),downbins,length(bwfiles[[2]])))
  
  minarray_downs_scaled_coll = apply(apply(minarray_downs, c(2,3), FUN = mean, na.rm = T),
                                     2, FUN = rev)
  
  
  #browser()
  #####stitch
  
  covmat_ups = (minarray_ups_scaled_coll + plusarray_ups_scaled_coll)/2
  covmat = (minarray_scaled_coll_corr + plusarray_scaled_coll)/2
  covmat_downs = (minarray_downs_scaled_coll + plusarray_downs_scaled_coll)/2
  
  final_matrix = t(rbind(covmat_ups, covmat, covmat_downs))
  rownames(final_matrix) = samnames
  
  return(final_matrix)
  
}




FIMO_parse_coords_stranded = function(regions, thresh, motif, 
                                      fimobfile = "fimo_bfile",
                                      countseq = F){
  
  r_plus = regions[strand(regions) == "+"]
  r_minus = regions[strand(regions) == "-"]
  
  hits_plus = runFimo(sequences = FASTA_w_coords(r_plus, write = F),
                      motifs = motif, 
                      parse_genomic_coord = T,
                      thresh = thresh, 
                      max_stored_scores = 100000,
                      norc = T,
                      bfile = fimobfile, 
                      meme_path = "/hpcnfs/home/ieo5559/meme/bin/")
  
  #this will need to be processed
  hits_minus = runFimo(sequences = FASTA_w_coords(r_minus, write = F, with_strand = F),
                       motifs = motif, 
                       parse_genomic_coord = F,
                       thresh = thresh, 
                       max_stored_scores = 100000,
                       norc = T,
                       bfile = fimobfile,
                       meme_path = "/hpcnfs/home/ieo5559/meme/bin/")
  
  
  ####first, take the end coordinate of the seqname of the hit, 
  ####then subtract the start coordinate of the hit. This is newend.
  ####new start is simply newend - width(hit)
  p1 = sapply(1:length(hits_minus), function(x){str_split(as.vector(seqnames(hits_minus))[x], "-")})
  
  newstart = as.numeric(sapply(p1, function(x){x[2]})) - end(hits_minus)
  newend = newstart + width(hits_minus)
  chrvec = sapply(str_split(sapply(p1, function(x){x[1]}), ":"), function(x){x[1]})
  
  hits_m_granges = GRanges(seqnames = chrvec, 
                           ranges = IRanges(start = newstart + 1, end = newend),
                           strand = rep("-", length(hits_minus))
  )
  
  if(countseq == T){
    print(length(unique(seqnames(runFimo(sequences = FASTA_w_coords(regions, write = F),
                                         motifs = motif, 
                                         parse_genomic_coord = F,
                                         thresh = thresh, 
                                         max_stored_scores = 100000,
                                         norc = T,
                                         bfile = fimobfile,
                                         meme_path = "/hpcnfs/home/ieo5559/meme/bin/")))))
  }
  
  mcols(hits_m_granges) = mcols(hits_minus)
  ##concatenate the two
  return(c(hits_plus, hits_m_granges))
  
}

FASTA_w_coords = function(regions, write = F, outfile, with_strand = F, useregnames = F){
  
  seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, regions)
  if(with_strand == T & useregnames == T){
    stop("both cannot be true")
  }
  if(with_strand == T){
    namesreg = paste0(seqnames(regions), ":", start(regions), "-", end(regions), "(", strand(regions), ")")
    names(seqs) = namesreg
    return(seqs)
  }else if(useregnames == T){
    names(seqs) = names(regions)
    return(seqs)
  }else{
    namesreg = paste0(seqnames(regions), ":", start(regions), "-", end(regions))
    names(seqs) = namesreg
    return(seqs)
  }
}


intex_ratio = function(bamfilelist, mane, pctx, genes_to_plot){
  
  #pass mane
  tx_filt = pctx[pctx$GENEID %in% genes_to_plot]
  
  ######counting
  
  n_exons = sapply(exonsBy(mane, by = "tx", use.names = T)[names(tx_filt)], length)
  tx_lengths = width(transcripts(mane, use.names = T)[names(tx_filt)])
  txdf = data.frame(row.names = names(tx_filt), "n_exons" = n_exons, "length" = tx_lengths)
  #txdf = txdf[txdf$length > 20000,]
  
  #longtx = txdf[txdf$length > 60000,]
  #exonsBy(mane, by = "tx", use.names = T)[rownames(longtx)]
  
  ## intron exon read ratio
  exons_all = exonsBy(mane, by = "tx", use.names = T)[rownames(txdf)]
  introns_all = intronsByTranscript(mane, use.names = T)[rownames(txdf)]
  
  # want multi-exonic
  exons = exons_all[sapply(introns_all, length) > 0]
  introns = introns_all[sapply(introns_all, length) > 0]
  
  exoncounts = summarizeOverlaps(features = exons, reads = bamfilelist, singleEnd=FALSE, 
                                 fragments=F, ignore.strand = F, BPPARAM = bpparam(),
                                 strandMode = 2)
  totexonlengths = sapply(sapply(exons, width), sum)
  
  introncounts = summarizeOverlaps(features = introns, reads = bamfilelist, singleEnd=FALSE, 
                                   fragments=F, ignore.strand = F, BPPARAM = bpparam(),
                                   strandMode = 2)
  totintronlengths = sapply(sapply(introns, width), sum)
  
  ######
  stopifnot(all.equal(rownames(exoncounts), 
            rownames(introncounts), 
            names(totexonlengths), 
            names(totintronlengths)))
  
  exon_rpks = sweep(assay(exoncounts), MARGIN = 1, FUN = "/",
                    STATS = totexonlengths)*1000
  intron_rpks = sweep(assay(introncounts), MARGIN = 1, FUN = "/",
                      STATS = totintronlengths)*1000
  
  return(list("exon_rpks" = exon_rpks, "intron_rpks" = intron_rpks))
  
}

spliceai_wrapper = function(inputseqs){
  
  if(is.null(names(inputseqs))){
    names(inputseqs) = 1:length(inputseqs)
  }
  
  writeXStringSet(inputseqs, "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/miscellany/splice/procap_windows.fa")
  
  system2(command = "python3", 
          args = c("/hpcnfs/data/GN2/gmandana/bin/spliceAIwrapper.py", 
                   "--fname /hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/miscellany/splice/procap_windows.fa", 
                   "--foutname /hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/miscellany/splice/output_donor_preds.txt"),
          env = )
  
  toreturn = read.table("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/miscellany/splice/output_donor_preds.txt")
  
  return(toreturn)
  
}

cov_without_collapse = function(bwfiles, regions, nbins, norm = F){
  
  if(is.null(names(regions))){
    names(regions) = paste0("id", 1:length(regions))
  }
  
  if(length(unique(width(regions))) != 1){
    stop("need to be same length")
  }
  
  samnames = basename(file_path_sans_ext(bwfiles[[1]]))
  
  r_p = regions[strand(regions) == "+"]
  r_m = regions[strand(regions) == "-"]
  
  plusmat = lapply(1:length(bwfiles[[1]]), function(x){
    bigwig_collect(bwfiles[[1]][x], regions = r_p, norm = F, nbins = nbins)
  })
  
  plusarray = array(unlist(plusmat), c(length(r_p),nbins,length(bwfiles[[1]])))
  
  ######now same for reverse
  minmat = lapply(1:length(bwfiles[[2]]), function(x){
    bigwig_collect(bwfiles[[2]][x], regions = r_m, norm = F, nbins = nbins)
  })
  
  minarray = array(unlist(minmat), c(length(r_m),nbins,length(bwfiles[[2]])))
  
  covmat = abind(plusarray, minarray, along = 1)
  
  return(covmat)
  
}


region_cleaning = function(regions, txdb, fdrthresh, 
                           marker_regions = NULL,
                           filterthreeprime = T,
                           lenfilt = NULL){
  
  ### this is for CSAW regions
  if(!is.null(marker_regions)){
    markers = tryCatch(unlist(marker_regions), error = function(e){marker_regions})
  }
  
  if(!is.null(lenfilt)){
    regions = regions[width(regions) > lenfilt]
  }
  
  filt = regions[regions$fdr >= fdrthresh & 
                   sapply(as.vector(seqnames(regions)), nchar) < 8]
  
  # use PROCAP, w/e
  if(!is.null(marker_regions)){
    filt = filt[countOverlaps(promoters(filt, upstream = 1000, downstream = 1000), 
                            slop(markers, 1000)) > 0]
  }
  
  if(filterthreeprime){
    filt = subsetByOverlaps(filt, unlist(threeUTRsByTranscript(txdb)), 
                            ignore.strand = F, invert = T, maxgap = 100)
  }
  
  return(filt)
  
}

nc_reg_annotation_MANE = function(mane, regs){
  
  txdb = mane
  
  if(is.null(names(regs))){
    names(regs) = paste0("id", 1:length(regs))
  }
  
  regs2 = promoters(regs, upstream = 0, downstream = 1)
  
  regs$annot = NA
  regs$class = NA
  
  maxd = 1000
  
  ##start pa-RNA
  olap1 = findOverlaps(invertStrand(promoters(regs2, upstream = 0, downstream = 1)), 
                       promoters(txdb, upstream = maxd*3, downstream = maxd*3), 
                       ignore.strand = F)
  
  annot1 = as.data.frame(olap1) %>% group_by(queryHits) %>% filter(row_number()==1)
  regs[annot1$queryHits]$annot = paste0(txdb[annot1$subjectHits]$V13, "-AS")
  regs[annot1$queryHits]$class = "pa-RNA"
  
  ## now rest
  unannot = regs[is.na(regs$annot)]
  
  olap2 = findOverlaps(resize(unannot, width = 10, fix = "end"), 
                       promoters(txdb, upstream = 1000, downstream = 0), 
                       ignore.strand = F)
  
  if(nrow(as.data.frame(olap2)) == 0){
    unannot2 = unannot
  }else{
    annot2 = as.data.frame(olap2) %>% group_by(queryHits) %>% filter(row_number()==1)
    
    unannot[annot2$queryHits]$annot = paste0(txdb[annot2$subjectHits]$V13, "-upstream_TSS")
    unannot[annot2$queryHits]$class = "upstreamTSS"
    unannot2 = unannot[is.na(unannot$annot)]
  }
  
  ### now termination defects
  
  olap3 = findOverlaps(resize(unannot2, width = 10, fix = "start"), 
                       resize(txdb, width = 10, fix = "end"), 
                       ignore.strand = F)
  
  if(nrow(as.data.frame(olap3)) == 0){
    unannot3 = unannot2
  }else{
    annot3 = as.data.frame(olap3) %>% group_by(queryHits) %>% filter(row_number()==1)
    
    unannot2[annot3$queryHits]$annot = paste0(txdb[annot3$subjectHits]$V13, "-termination_defect")
    unannot2[annot3$queryHits]$class = "termination_defect"
    unannot3 = unannot2[is.na(unannot2$annot)]
  }
  
  ## now intronic enhancer
  olap4 = findOverlaps(invertStrand(promoters(unannot3, upstream = 0, downstream = 1)), 
                       txdb, 
                       ignore.strand = F)
  
  if(nrow(as.data.frame(olap4)) == 0){
    unannot4 = unannot3
  }else{
    annot4 = as.data.frame(olap4) %>% group_by(queryHits) %>% filter(row_number()==1)
    
    unannot3[annot4$queryHits]$annot = paste0(txdb[annot4$subjectHits]$V13, "-intronic_enhancer")
    unannot3[annot4$queryHits]$class = "intronic_enhancer"
    unannot4 = unannot3[is.na(unannot3$annot)]
  }
  
  ###the rest are probs enhancer
  ###find the ones that are tails of the other
  
  #intermed = c(unannot3[!is.na(unannot3$annot)],
  #             unannot2[!is.na(unannot2$annot)],
  #             unannot[!is.na(unannot$annot)],
  #             regs[!is.na(regs$annot)])
  #
  #olap6 = findOverlaps(promoters(unannot4, upstream = 0, downstream = 1),
  #                     promoters(resize(intermed, width = 1, fix = "end"),
  #                               upstream = 0, downstream = 5000), ignore.strand = F)
  #
  #annot6 = as.data.frame(olap6) %>% group_by(queryHits) %>% filter(row_number()==1)
  #
  #unannot5[annot6$queryHits]$annot = paste0(intermed[annot6$subjectHits]$annot, "_2")
  #unannot5[annot6$queryHits]$class = intermed[annot6$subjectHits]$class
  
  ###rest are enhancer
  
  unannot4[is.na(unannot4$annot)]$class = "enhancer"
  unannot4[is.na(unannot4$annot)]$annot = "enhancer"
  
  #combine them all
  
  all = c(unannot4[!is.na(unannot4$annot)],
          unannot3[!is.na(unannot3$annot)],
          unannot2[!is.na(unannot2$annot)],
          unannot[!is.na(unannot$annot)],
          regs[!is.na(regs$annot)])
  
  siW_annotated_2 = all[names(regs)]
  
  return(siW_annotated_2)
  
}

BL = import("/hpcnfs/data/GN2/gmandana/annotation/small_RNA_tothrow_and_just_tRNA_hg38_encode_BL_bed4.bed")

####porco Dio
## I need a function to generate a metaplot...
## given the coverage matrix, and metadata  
metaplotter = function(matrix, metadata, return_df = F, custom_x_ticks = NULL){
  
  ##metadata contains replicate info, i.e., which cols to calculate mean and SD over
  #### we need to do some averaging
  avging_mat = as.data.frame(cbind(1:length(metadata), sapply(str_split(metadata, "_R[0-9]"), function(x){x[[1]]})))
  
  rows_to_avg = split(1:length(metadata), avging_mat$V2)
  
  # the trycatch is for samples with only one replicate
  avg_prof = sapply(rows_to_avg, function(x){
    tryCatch(colMeans(matrix[x,], na.rm = T), error = function(e){matrix[x,]})
    })
  
  #sd_vals = sapply(rows_to_avg, function(x){
  #  tryCatch(colSds(matrix[x,], na.rm = T), error = function(e){rep(0, dim(matrix)[2])})
  #})
  ##SEM
  sd_vals = sapply(rows_to_avg, function(x){
    tryCatch(colSds(matrix[x,], na.rm = T)/sqrt(length(x)), error = function(e){rep(0, dim(matrix)[2])})
  })
  
  melted_prof = melt(avg_prof); colnames(melted_prof) = c('bin', 'group', 'value')
  melted_sds = melt(sd_vals); colnames(melted_sds) = c('bin_sd', 'group_sd', 'value_sd')
  final = cbind(melted_prof, melted_sds)
  
  if(!return_df){
    metaplot = final %>% 
      ggplot(aes(x = bin, y = value, fill = group)) + 
      geom_line(aes(color = group)) + 
      geom_ribbon(aes(y = value, ymin = value+(value_sd), 
                      ymax = value-(value_sd), group = group), 
                  alpha = 0.5) + 
      theme_classic() + xlab("") + ylab("")
    
    return(metaplot)
  }else{
    return(final)
  }
  
}

#see https://github.com/zanglab/SICER2/blob/master/sicer/src/compare_two_libraries_on_islands.py
BH_FDR_calc = function(bedfile){
  
  pvals = bedfile$V7
  FDR = (pvals*length(pvals))/rank(pvals)
  
  FDR[FDR > 1] = 1
  FDR[is.na(FDR)] = 1;
  return(FDR)
  
}

SICER_region_cleaning = function(bedfile, fdrthresh, mane, 
                                 gencode, k27, PROCAP, lenfilt){
  
  bedfile$FDR = BH_FDR_calc(bedfile)
  
  filt1 = bedfile[bedfile$FDR < fdrthresh]
  filt2 = filt1[names(subsetByOverlaps(promoters(filt1, upstream = 1000, downstream = 0),
                                       unlist(threeUTRsByTranscript(mane)), 
                                       ignore.strand = F, invert = T, maxgap = 0))]
  
  ### you are asked to remove SNHG
  SNHG = gencode[grepl("SNHG", gencode$V13)]
  
  filt3 = subsetByOverlaps(filt2, SNHG, invert = T)
  
  # just to double check
  filt3$k27 = countOverlaps(promoters(filt3, upstream = 1000, downstream = 1000), slop(PROCAP, 1000)) > 0
  filt3$procap = countOverlaps(promoters(filt3, upstream = 1000, downstream = 1000), slop(k27, 1000)) > 0
  
  filt4 = filt3[filt3$procap == T | filt3$k27 == T]
  filt5 = filt4[width(filt4) >= 2000]
  
  filt5$procap = filt5$k27 = NULL
  
  return(filt5[order(filt5$FDR, decreasing = F)])
}

SICER_region_splitting = function(regs, bamfiles, PROCAP, 
                                  pvalthresh, selTSS = NULL){
  
  if(is.null(names(regs))){
    stop("regions must be named")
  }
  if(is.null(names(PROCAP))){
    stop("PROCAP must be named")
  }
  
  if(is.null(selTSS)){
    rel_procap2 = subsetByOverlaps(subsetByOverlaps(PROCAP, regs),
                                   promoters(regs, upstream = 0, downstream = 5000),
                                   invert = T)
    rel_procap_ups2 = summarizeOverlaps(features = promoters(rel_procap2, 
                                                             upstream = 1000,
                                                             downstream = 0), 
                                        reads = bamfiles, 
                                        inter.feature = FALSE,
                                        singleEnd=FALSE, 
                                        fragments=F, ignore.strand = F, 
                                        BPPARAM = MulticoreParam(),
                                        strandMode = 2)
    rel_procap_downs2 = summarizeOverlaps(features = promoters(rel_procap2, 
                                                               upstream = 0,
                                                               downstream = 1000), 
                                          reads = bamfiles, 
                                          inter.feature = FALSE,
                                          singleEnd=FALSE, 
                                          fragments=F, ignore.strand = F, 
                                          BPPARAM = MulticoreParam(),
                                          strandMode = 2)
    
    counts_PROCAP2 = cbind(assay(rel_procap_ups2), assay(rel_procap_downs2))
    colnames(counts_PROCAP2) = c("up_r1", "up_r2", "down_r1", "down_r2")
    
    desmat = matrix(rep(c(1,0,1,0,0,1,0,1), 2),4,2, byrow = T)
    colnames(desmat) = c("up","down")
    
    dlist = asDGEList(SummarizedExperiment(counts_PROCAP2), assay.id=1)
    dlist = estimateDisp(dlist, desmat)
    fit = glmQLFit(dlist, desmat, robust=T)
    results = glmQLFTest(fit, contrast = makeContrasts(down-up, levels = desmat))
    results$table$pval = -log10(results$table$PValue)
    
    thresh = pvalthresh
    sel_tss = rel_procap2[rownames(results$table[results$table$pval > thresh & 
                                                   results$table$logFC > 0,])]
  }else{
    sel_tss = selTSS
  }
  
  
  ########
  ######## use sel TSS to split
  sel_tss_point = resize(sel_tss, width = 1, fix = "center")
  sel_tss_point_2 = resize(reduce(sel_tss_point, min.gapwidth = 1000), width = 1, fix = "center")
  to_correct = subsetByOverlaps(regs, sel_tss_point_2)
  
  split_up = lapply(1:length(to_correct), function(x){
    break_tx(to_correct[x], sel_tss_point_2)
  })
  
  unmerged_nctx2 = c(regs[!(names(regs) %in% names(to_correct))], 
                     do.call('c', split_up))
  
  return(unmerged_nctx2)
  
}

get_first_SS_info = function(ranges, juncs){
  
  ol = as.data.frame(findOverlaps(juncs, ranges))
  ol$strand = as.vector(strand(ranges[ol$subjectHits]))
  firstSSinds = ol %>%
    group_by(subjectHits) %>%
    mutate(isfirst = case_when(strand == "+" ~ row_number() == 1,
                               .default = row_number() == n())) %>%
    dplyr::filter(isfirst == T)
  
  firstSS = juncs[firstSSinds$queryHits]
  
  ######### DONOR SITES
  #dist
  dist_5 = distance(resize(firstSS,width = 1, fix = "start"), 
           resize(ranges[firstSSinds$subjectHits], width = 1, fix = "start"))
  
  #maxent
  first_donor_ss_win_9 = shiftStranded(resize(resize(firstSS, fix = "start", width = 1), 
                                              fix = "center", width = 9), 1)
  first_donor_ss_win_9_seqs = suppressWarnings(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, first_donor_ss_win_9))
  
  Maxent_donor = sapply(1:length(first_donor_ss_win_9_seqs), function(x){
    calculateMaxEntScanScore(first_donor_ss_win_9_seqs[x], 5)
  })
  
  Maxent_donor_scores = as.numeric(Maxent_donor)
  
  ####### ACCEPTOR SITES
  threeP = resize(juncs, width = 1, fix = "end")
  threeP_ordered = makeGRangesFromDataFrame(as.data.frame(threeP) %>% 
                                              group_by(seqnames) %>% 
                                              arrange(start, .by_group = T), 
                                            keep.extra.columns = T)
  
  ol_3p = as.data.frame(findOverlaps(threeP_ordered, ranges))
  ol_3p$strand = as.vector(strand(ranges[ol_3p$subjectHits]))
  first_acc_SSinds = ol_3p %>%
    group_by(subjectHits) %>%
    mutate(isfirst = case_when(strand == "+" ~ row_number() == 1,
                               .default = row_number() == n())) %>%
    dplyr::filter(isfirst == T)
  
  first_accpetor = threeP_ordered[first_acc_SSinds$queryHits]
  
  #dist
  dist_3 = distance(first_accpetor, resize(ranges[first_acc_SSinds$subjectHits], width = 1, fix = "start"))
  
  #maxent
  first_accpetor_win = promoters(first_accpetor, upstream = 19, downstream = 4)
  first_accpetor_win_seqs = suppressWarnings(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, first_accpetor_win))
  
  Maxent_acceptor = sapply(1:length(first_accpetor_win_seqs), function(x){
    calculateMaxEntScanScore(first_accpetor_win_seqs[x], 3)
  })
  
  Maxent_acceptor_scores = as.numeric(Maxent_acceptor)
  
  return_obj = list("firstSS" = firstSS,
                    "reltx" = ranges[firstSSinds$subjectHits],
                    "first_donor_seqs" = first_donor_ss_win_9_seqs,
                    "MaxEnt_donor" = Maxent_donor_scores,
                    "first_accpetor" = first_accpetor,
                    "first_acceptor_seqs" = first_accpetor_win_seqs,
                    "Maxent_acceptor" = Maxent_acceptor_scores,
                    "dist_to_5" = dist_5,
                    "dist_to_3" = dist_3)
  
  return(return_obj)
  
}






















