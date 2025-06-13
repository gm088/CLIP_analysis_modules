library(dplyr)
library(tidyr)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(GenomicRanges)
library(slider)
library(csaw)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(universalmotif)
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")
library(kmer)
library(reshape2)
library(pheatmap)
library(dendextend)
library(kernlab)
library(e1071)
library(factoextra)
library(FactoMineR)
library(gtools)
library(memes)
library(caret)
library(rminer)

####spaced kmer counts

generate_occurence_matrix = function(seqs, max.spacer = 10, kmerlength = 3, hits_per_seq = T){
  
  options(expressions=1e5)  #needed for gtools::permutations for n > 40
  
  #generate all possible kmers from the DNA alphabet
  all_kmers = apply(permutations(4, kmerlength, v = c("A", "C", "G", "T"), repeats.allowed = T), 
                        1, paste, collapse = "")
  #generate all possible permutations of kmers
  permuts = permutations(64, 2, all_kmers, repeats.allowed = T)
  permuts_tags = sapply(1:nrow(permuts), function(x){paste(permuts[x,], collapse = '..')})
  
  if(hits_per_seq == T){
    #count OCCURENCE ONLY of a kmer-kmer combination, spacing between 0 and max.spacer
    occurence_matrix = do.call('cbind', lapply(0:max.spacer, function(y){
      sapply(1:nrow(permuts), function(x){
        #construct regex for spaced kmers
        sum(sapply(str_locate_all(seqs, paste(permuts[x,], collapse = sprintf('.{%d}', y))), nrow))
      })
    })
    )
  }else if(hits_per_seq == F){
    #count number of seqs with kmer-kmer combination, spacing between 0 and max.spacer
    occurence_matrix = do.call('cbind', lapply(0:max.spacer, function(y){
      sapply(1:nrow(permuts), function(x){
        #construct regex for spaced kmers
        sum(as.numeric(sapply(
          str_locate_all(seqs, paste(permuts[x,], collapse = sprintf('.{%d}', y))), nrow) > 0))
      })
    })
    )
  }
  
  rownames(occurence_matrix) = permuts_tags
  colnames(occurence_matrix) = 0:max.spacer
  return(occurence_matrix)
}

spaced_kmer = function(regions, max.spacer = 10, kmerlength = 3){
  
  options(expressions=1e5)  #needed for gtools::permutations for n > 40
  
  seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, regions, as.character = T)
  seqs = gsub("N", "A", seqs) #I did this only becuase there is ONE N in ONE sequence
  
  #generate all possible kmers from the DNA alphabet
  top_kmers = rownames(as.data.frame(colSums(kcount(str_split(seqs, ""),kmerlength))) %>% 
                         arrange(desc(colSums(kcount(str_split(seqs, ""), kmerlength)))) %>% 
                         head(20))
  #generate all possible permutations of kmers
  permuts = permutations(20, 2, top_kmers, repeats.allowed = T)
  #permuts_tags = sapply(1:nrow(permuts), function(x){paste(permuts[x,], collapse = '..')})
  
  allregex = unlist(lapply(1:nrow(permuts), function(x){
    sapply(1:max.spacer, function(y){paste(permuts[x,], collapse = sprintf('.{%d}', y))})
    }))
  
  #count occurence of kmer-kmer combination in each sequence, spacing between 0 and max.spacer
  occurence_matrix = sapply(allregex, function(x){
    #construct regex for spaced kmers
    sapply(str_locate_all(seqs, x), nrow)
    })
  
  rownames(occurence_matrix) = names(regions)
  return(occurence_matrix)
}

customimportbed6 = function(fname){
  
  regions = makeGRangesFromDataFrame(read.table(fname, header = F), 
                           keep.extra.columns = T, ignore.strand = F,
                           seqnames.field = "V1",start.field = "V2",
                           end.field = "V3",strand.field = "V6")
  
  return(regions)
}

homogenise_region = function(ranges, width_to_make = 5000, id_prefix = ""){
  
  peaks_hom = resize(ranges, width = width_to_make, fix = "center")
  peaks_hom = peaks_hom[!duplicated(csaw::mergeWindows(peaks_hom, tol = 0)$ids)]
  names(peaks_hom) = paste(id_prefix,"_id", 1:length(peaks_hom), sep = "")
  
  return(peaks_hom)
}

custom_homogenise = function(regions, ups, downs, id_prefix = ""){
  
  extended = promoters(regions, upstream = ups, downstream = downs)
  extended = extended[!duplicated(csaw::mergeWindows(extended, tol = 0)$ids)]
  names(extended) = paste(id_prefix,"_id", 1:length(extended), sep = "")
  return(extended)
}

kmercount = function(k, regions){
  
  ##wrapper around kcount because it needs that annoying string split step
  seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, regions, as.character = T)
  #seqs = seqs[grep("N", seqs, invert = T)]
  seqs = gsub("N", "A", seqs) #I did this only becuase there is ONE N in ONE sequence
  seqs_rawstr = sapply(seqs, str_split, "")
  kmers = kcount(seqs_rawstr, k)
  
  return(kmers)
}

process_fimo_output = function(ranges, fimo_granges, sum = T){
  
  if(sum == T){
    prepr = as.data.frame(fimo_granges) %>% 
      dplyr::select(-c(motif_id, score, qvalue)) %>% ##remove useless columns
      mutate(score = -log10(pvalue)) %>% 
      group_by(seqnames) %>% 
      group_by(motif_alt_id, .add = T) %>% 
      summarise(sumscore = sum(score)) %>% ##take the sum of all pvalues of hits
      pivot_wider(names_from = seqnames, values_from = c(sumscore))
  }else{
    prepr = as.data.frame(fimo_granges) %>% 
      dplyr::select(-c(motif_id, score, qvalue)) %>% 
      mutate(score = -log10(pvalue)) %>% 
      group_by(seqnames) %>% 
      group_by(motif_alt_id, .add = T) %>% 
      summarise(maxscore = max(score)) %>% ## take only the best hit
      pivot_wider(names_from = seqnames, values_from = c(maxscore))
  }
  processed = t(as.data.frame(prepr)[,2:ncol(prepr)]) #because the first column becomes
                                                        #motif names
  colnames(processed) = prepr$motif_alt_id
  processed = as.data.frame(processed)
  
  ###in case there are some sequences without hits
  if(length(ranges) != length(rownames(processed))){
    missing = names(ranges)[!(names(ranges) %in% rownames(processed))]
    toadd = data.frame(row.names = missing, 
                       matrix(data = NA, nrow = length(missing), ncol = ncol(processed), byrow = FALSE, 
                              dimnames = NULL))
    colnames(toadd) = prepr$motif_alt_id
    
    final = rbind(processed, toadd)
    
    final[is.na(final)] = 0
    return(as.data.frame(final))
  }else{
    processed[is.na(processed)] = 0
    return(processed)
  }
  
}

repeatcounts = function(ranges, repeatsbed){
  
  #I want a dataframe with seqs, and num hits
  families = unique(repeatsbed$V4)
  
  df = do.call('cbind', lapply(families, function(x){
    countOverlaps(ranges, repeatsbed[repeatsbed$V4 == x])
    }))
  colnames(df) = families
  
  return(df[,colSums(df) > 0])
  
}

GC_fraction = function(ranges){
  
  GC_frac = letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, ranges), 
                            letters = c("AT", "GC"), as.prob = T)[,2]
  GC_df = as.data.frame(GC_frac)
  rownames(GC_df) = names(ranges)
  return(GC_df)
  
}

region_tiling = function(regions, id_prefix = ""){
  
  window = promoters(regions, upstream = 1000, downstream = 2000)
  tosplit = tile(window, 6)
  
  #I could not find a better way to index into a GRangesList
  #for the minus strand, have to extract last element first
  win1 = do.call('c', lapply(1:length(window), function(x){if(as.vector(strand(window[x])) == "+"){
    unlist(tosplit[x,])[1]}else{unlist(tosplit[x,])[6]}}))
  win2 = do.call('c', lapply(1:length(window), function(x){if(as.vector(strand(window[x])) == "+"){
    unlist(tosplit[x,])[2]}else{unlist(tosplit[x,])[5]}}))
  win3 = do.call('c', lapply(1:length(window), function(x){if(as.vector(strand(window[x])) == "+"){
    unlist(tosplit[x,])[3]}else{unlist(tosplit[x,])[4]}}))
  win4 = do.call('c', lapply(1:length(window), function(x){if(as.vector(strand(window[x])) == "+"){
    unlist(tosplit[x,])[4]}else{unlist(tosplit[x,])[3]}}))
  win5 = do.call('c', lapply(1:length(window), function(x){if(as.vector(strand(window[x])) == "+"){
    unlist(tosplit[x,])[5]}else{unlist(tosplit[x,])[2]}}))
  win6 = do.call('c', lapply(1:length(window), function(x){if(as.vector(strand(window[x])) == "+"){
    unlist(tosplit[x,])[6]}else{unlist(tosplit[x,])[1]}}))
  
  all = list(win1, win2, win3, win4, win5, win6)
  
  #remember to name them
  for(i in 1:6){
    names(all[[i]]) = paste0(id_prefix, "_id", 1:length(regions))
  }
  return(all)
  
}

bigwig_collect = function(bwfile, regions, label, nbins = NULL){
  
  prof = do.call('rbind', coverageFromBW(bwfile = bwfile, gr = regions, Nbin = nbins))
  
  ##reverse those that have minus strand
  ###some bullshit fix
  prof[as.vector(strand(regions)) == "-",] = t(apply(prof[as.vector(strand(regions)) == "-",], 1, rev))
  
  colnames(prof) = paste0(label, "_", 1:nbins)
  return(nainorm(prof))
  
}

collect_and_train = function(fg_regions, bg_regions, identifier, motifs, repeats){
  
  hits_fg = runFimo(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, fg_regions), motifs = motifs, norc = T, thresh = 1e-4)
  hits_bg = runFimo(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, bg_regions), motifs = motifs, norc = T, thresh = 1e-4)

  fimo_fg = process_fimo_output(fg_regions, hits_fg)
  fimo_bg = process_fimo_output(bg_regions, hits_bg)

  ### kmer counts per base position
  
  fg_threemers = kmercount(3, fg_regions)
  bg_threemers = kmercount(3, bg_regions)

  fg_fourmers = kmercount(4, fg_regions)
  bg_fourmers = kmercount(4, bg_regions)

  #fg_sixmers = kmercount(6, fg_hom)
  #bg_sixmers = kmercount(6, bg_hom)

  ## GC content
  
  fg_GC = GC_fraction(fg_regions)
  bg_GC = GC_fraction(bg_regions)

  ### repeats
  
  fg_repeats = repeatcounts(fg_regions, repeats)
  bg_repeats = repeatcounts(bg_regions, repeats)

  ### collect all the features
  fg_all = bind_cols(fimo_fg, as.data.frame(fg_threemers), as.data.frame(fg_fourmers), fg_GC, fg_repeats)
  bg_all = bind_cols(fimo_bg, as.data.frame(bg_threemers), as.data.frame(bg_fourmers), bg_GC, bg_repeats)

  
  ###train SVM
  
  combined = bind_rows(fg_all, bg_all)
  combined$group = as.factor(
    sapply(str_split(rownames(combined), "_"), function(x){x[1]})
  )
  combined[is.na(combined)] = 0
  
  smp = bind_rows(fg_all[sample(1:nrow(fg_all), 300, replace = F),], 
                  bg_all[sample(1:nrow(bg_all), 300, replace = F),])
  
  smp$group = as.factor(
    sapply(str_split(rownames(smp), "_"), function(x){x[1]})
  )
  smp[is.na(smp)] = 0
  
  gammas = 2^(-15:-1)
  costs = 2^(-10:10)
  
  svmgs <- tune(svm,
                train.x = smp[-ncol(smp)],
                train.y = smp[,ncol(smp)],
                type = "C-classification",
                kernel = "radial", 
                scale = T,
                ranges = list(gamma = gammas, cost = costs),
                tunecontrol = tune.control(cross = 4)
  )
  
  saveRDS(svmgs, file = paste0("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/files_for_svm/model_",
                               identifier, ".RDS"))
  
  result <- as.data.frame(predict(svmgs$best.model, combined[-ncol(combined)]))
  result$group = combined[,ncol(combined)]
  
  confusionMatrix(result$`predict(svmgs$best.model, combined[-ncol(combined)])`, 
                  result$group, positive = "fg")
  
  
  svm_imp = Importance(svmgs$best.model, 
                       data = smp[-ncol(smp)],
                       PRED = predict)
  
  
}

svm_train = function(fg, bg, train_set_size, kernel, gammas, costs){
  
  ###train SVM
  combined = bind_rows(fg, bg)
  combined$group = as.factor(
    sapply(str_split(rownames(combined), "_"), function(x){x[1]})
  )
  combined[is.na(combined)] = 0
  
  smp = bind_rows(fg[sample(1:nrow(fg), train_set_size, replace = F),], 
                  bg[sample(1:nrow(bg), train_set_size, replace = F),])
  
  smp$group = as.factor(
    sapply(str_split(rownames(smp), "_"), function(x){x[1]})
  )
  smp[is.na(smp)] = 0
  
  gammas = gammas
  costs = costs
  
  svmgs <- tune(svm,
                train.x = smp[-ncol(smp)],
                train.y = smp[,ncol(smp)],
                type = "C-classification",
                kernel = kernel, 
                scale = T,
                ranges = list(gamma = gammas, cost = costs),
                tunecontrol = tune.control(cross = 5)
  )
  test_set = combined[!(rownames(combined) %in% rownames(smp)),]
  
  result <- as.data.frame(predict(svmgs$best.model, test_set[-ncol(test_set)]))
  result$group = test_set[,ncol(test_set)]
  
  consmat = confusionMatrix(result$`predict(svmgs$best.model, test_set[-ncol(test_set)])`, 
                  result$group, positive = "fg")
  
  
  svm_imp = Importance(svmgs$best.model, 
                       data = smp[-ncol(smp)],
                       PRED = predict)
  
  return(list(svmgs, svm_imp, consmat, smp, test_set))
  
}

collect = function(fg, bg, motifs, repeats, repeats2, expand = T,
                   upstream = 500, downstream = 1000, fimo_bg_file, fimo_thresh,
                   fimo_max_rep, histonebw1,
                   histonebw2){

  if(expand == T){
    fg_regions = homogenise_region(fg, downstream, "fg")
    bg_regions = homogenise_region(bg, downstream, "bg")
  }else{
    fg_regions = custom_homogenise(fg, upstream, downstream, "fg")
    bg_regions = custom_homogenise(bg, upstream, downstream, "bg")
  }
  
  hits_fg = runFimo(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, fg_regions), 
                    motifs = motifs, norc = T, thresh = fimo_thresh, bfile = fimo_bg_file,
                    max_stored_scores = fimo_max_rep)
  hits_bg = runFimo(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, bg_regions), 
                    motifs = motifs, norc = T, thresh = fimo_thresh, bfile = fimo_bg_file,
                    max_stored_scores = fimo_max_rep)
  
  fimo_fg = process_fimo_output(fg_regions, hits_fg, sum = F)
  fimo_bg = process_fimo_output(bg_regions, hits_bg, sum = F)
  
  ##histone mark
  #fg_hist1 = bigwig_collect(bwfile = histonebw1[[1]], regions = fg_regions, 
  #                          histonebw1[[2]], nbins = 30)
  #bg_hist1 = bigwig_collect(bwfile = histonebw1[[1]], regions = bg_regions, 
  #                          histonebw1[[2]], nbins = 30)
  #
  #fg_hist2 = bigwig_collect(bwfile = histonebw2[[1]], regions = fg_regions, 
  #                          histonebw2[[2]], nbins = 30)
  #bg_hist2 = bigwig_collect(bwfile = histonebw2[[1]], regions = bg_regions, 
  #                          histonebw2[[2]], nbins = 30)
  
  ### kmer counts per base position
  
  fg_threemers = kmercount(3, fg_regions)
  bg_threemers = kmercount(3, bg_regions)
  
  fg_fourmers = kmercount(4, fg_regions)
  bg_fourmers = kmercount(4, bg_regions)
  
  fg_sixmers = kmercount(6, fg_regions)
  bg_sixmers = kmercount(6, bg_regions)
  
  ## GC content
  fg_GC = GC_fraction(fg_regions)
  bg_GC = GC_fraction(bg_regions)
  
  ### repeats
  #fg_repeats = repeatcounts(fg_regions, repeats)
  #bg_repeats = repeatcounts(bg_regions, repeats)
  
  #fg_simple_repeats = repeatcounts(fg_regions, repeats2)
  #bg_simple_repeats = repeatcounts(bg_regions, repeats2)
  
  ### collect all the features
  fg_all = bind_cols(fimo_fg, 
                     as.data.frame(fg_threemers), 
                     as.data.frame(fg_fourmers),
                     as.data.frame(fg_sixmers),
                     fg_GC)
  
  bg_all = bind_cols(fimo_bg, 
                     as.data.frame(bg_threemers), 
                     as.data.frame(bg_fourmers), 
                     as.data.frame(bg_sixmers),
                     bg_GC)

  return(list(fg_all, bg_all))
  
}
  




