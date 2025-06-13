library(fitdistrplus)

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






