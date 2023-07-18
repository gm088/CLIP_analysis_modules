library(fitdistrplus)
library(AnnotationDbi)

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














