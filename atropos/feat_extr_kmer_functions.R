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
library(fmsb)
library(ROCit)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

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
  extended = extended[!duplicated(csaw::mergeWindows(extended, tol = 0, ignore.strand=FALSE)$ids),]
  names(extended) = paste(id_prefix,"_id", 1:length(extended), sep = "")
  return(extended)
}

stratify_by_score = function(regions, col_to_strat, num.quants = 4){

  labels = 1:num.quants
  regions$quantile = cut(as.data.frame(regions)[,col_to_strat], 
                            breaks = rev(quantile(as.data.frame(regions)[,col_to_strat], 
                                                  probs = seq(0,1,1/num.quants))), 
                            labels = labels, 
                            include.lowest = T)
  
  regions$quantile = factor(regions$quantile, levels = c(5,4,3,2,1))
  
  list = lapply(rev(labels), function(x){regions[regions$quantile == x]})
  return(list)
}

fimo_bg_generator = function(ranges, outfile){
  
  #FIMO uses a zero-order model
  freqs = letterFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, ranges), 
                          letters = c("A", "C", "G", "T"), as.prob = T)
  
  write("#zero order", outfile)
  write(paste("A", colMeans(freqs)[1]), outfile, append = T, ncolumns = 1000)
  write(paste("C", colMeans(freqs)[2]), outfile, append = T, ncolumns = 1000)
  write(paste("G", colMeans(freqs)[3]), outfile, append = T, ncolumns = 1000)
  write(paste("T", colMeans(freqs)[4]), outfile, append = T, ncolumns = 1000)
  
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

bigwig_collect = function(bwfile, regions, label, norm = F, nbins = NULL){
  
  prof = do.call('rbind', coverageFromBW(bwfile = bwfile, gr = regions, Nbin = nbins))
  
  ##reverse those that have minus strand
  ###some bullshit fix
  prof[as.vector(strand(regions)) == "-",] = t(apply(prof[as.vector(strand(regions)) == "-",], 1, rev))
  
  #colnames(prof) = paste0(label, "_", 1:nbins)
  if(norm){
    return(nainorm(prof))
  }else{
    return(prof)
  }
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

svm_train = function(fg, bg, train_set_size, kernel, gammas, costs, w_prob = TRUE){
  
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
                probability = w_prob,
                ranges = list(gamma = gammas, cost = costs),
                tunecontrol = tune.control(cross = 5)
  )
  test_set = combined[!(rownames(combined) %in% rownames(smp)),]
  
  result <- as.data.frame(predict(svmgs$best.model, test_set[-ncol(test_set)]))
  result$group = test_set[,ncol(test_set)]
  probs = predict(svmgs$best.model, test_set[-ncol(test_set)], probability = T)
  
  consmat = confusionMatrix(result$`predict(svmgs$best.model, test_set[-ncol(test_set)])`, 
                  result$group, positive = "fg")
  
  
  svm_imp = Importance(svmgs$best.model, 
                       data = smp[-ncol(smp)],
                       PRED = predict)
  
  ROC = rocit(attr(probs, "probabilities")[,1], test_set$group)
  
  return(list("model" = svmgs, "imp_est" = svm_imp, 
              "constmat" = consmat, "train_set" = smp, "val_set" = test_set, 
              "predictions" = result, "probs" = probs, "roc_obj" = ROC))
  
}

write_tsv_output_for_python = function(fg_mat, bg_mat, outfile){
  
  combined = bind_rows(fg_mat, bg_mat)
  combined$group = as.factor(
    sapply(str_split(rownames(combined), "_"), function(x){x[1]})
  )
  combined[is.na(combined)] = 0
  
  write.table(combined, file = outfile,
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = T)
  
}

collect = function(fg, bg, RBP_motifs, TF_motifs, repeats, repeats2, expand = T,
                   upstream = 500, downstream = 1000, fimo_thresh,
                   fimo_max_rep, sixmer = F, alreadydefinedregions = F,
                   FIMO_TF_norc = F){

  if(alreadydefinedregions == F){
    if(expand == T){
      print(paste0("omogenise region: "), downstream)
      fg_regions = homogenise_region(fg, downstream, "fg")
      bg_regions = homogenise_region(bg, downstream, "bg")
    }else{
      print(paste0("custom homogenise region: "), upstream, " ", downstream)
      fg_regions = custom_homogenise(fg, upstream, downstream, "fg")
      bg_regions = custom_homogenise(bg, upstream, downstream, "bg")
    }
  }else{
    print("using provided regions")
    fg_regions = fg
    bg_regions = bg
  }
  
  #make fimo bg for fg and bg
  print("making FIMO bg files...")
  fimo_bg_generator(fg_regions, "fg_order0.param")
  fimo_bg_generator(bg_regions, "bg_order0.param")
  
  ###RBP motifs
  print("starting FIMO_RBP...")
  RBP_hits_fg = runFimo(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, fg_regions), 
                    motifs = RBP_motifs, norc = T, thresh = fimo_thresh, bfile = "fg_order0.param",
                    max_stored_scores = fimo_max_rep)
  
  RBP_hits_bg = runFimo(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, bg_regions), 
                    motifs = RBP_motifs, norc = T, thresh = fimo_thresh, bfile = "bg_order0.param",
                    max_stored_scores = fimo_max_rep)
  
  RBP_fimo_fg = process_fimo_output(fg_regions, RBP_hits_fg, sum = F)
  RBP_fimo_bg = process_fimo_output(bg_regions, RBP_hits_bg, sum = F)
  
  
  ###plot output of this for checking
  pdf("fimo_output_summary.pdf")
  par(mfrow = c(2,2))
  
  hist(-log10(RBP_hits_fg$pvalue), breaks = 100)
  hist(melt(RBP_fimo_fg)$value, breaks = 100)
  
  hist(-log10(RBP_hits_bg$pvalue), breaks = 100)
  hist(melt(RBP_fimo_bg)$value, breaks = 100)
  
  ###TF motifs
  print("starting FIMO_TF...")
  TF_hits_fg = runFimo(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, fg_regions), 
                        motifs = TF_motifs, norc = FIMO_TF_norc, thresh = fimo_thresh, bfile = "fg_order0.param",
                        max_stored_scores = fimo_max_rep)
  
  TF_hits_bg = runFimo(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, bg_regions), 
                        motifs = TF_motifs, norc = FIMO_TF_norc, thresh = fimo_thresh, bfile = "bg_order0.param",
                        max_stored_scores = fimo_max_rep)
  
  TF_fimo_fg = process_fimo_output(fg_regions, TF_hits_fg, sum = F)
  TF_fimo_bg = process_fimo_output(bg_regions, TF_hits_bg, sum = F)
  
  par(mfrow = c(2,2))
  
  hist(-log10(TF_hits_fg$pvalue), breaks = 100)
  hist(melt(TF_fimo_fg)$value, breaks = 100)
  
  hist(-log10(TF_hits_bg$pvalue), breaks = 100)
  hist(melt(TF_fimo_bg)$value, breaks = 100)
  
  #stop plotting
  dev.off()
  
  ## GC content
  fg_GC = GC_fraction(fg_regions)
  bg_GC = GC_fraction(bg_regions)
  
  ### kmer counts per base position
  print("counting kmers...")
  fg_threemers = kmercount(3, fg_regions)
  bg_threemers = kmercount(3, bg_regions)
  
  fg_fourmers = kmercount(4, fg_regions)
  bg_fourmers = kmercount(4, bg_regions)
  
  if(sixmer == T){
    
    print("counting 6mers...")
    fg_sixmers = kmercount(6, fg_regions)
    bg_sixmers = kmercount(6, bg_regions)
    
    ### collect all the features
    fg_all = bind_cols(RBP_fimo_fg, 
                       TF_fimo_fg,
                       as.data.frame(fg_threemers), 
                       as.data.frame(fg_fourmers),
                       as.data.frame(fg_sixmers),
                       fg_GC)
    
    bg_all = bind_cols(RBP_fimo_bg, 
                       TF_fimo_bg, 
                       as.data.frame(bg_threemers), 
                       as.data.frame(bg_fourmers), 
                       as.data.frame(bg_sixmers),
                       bg_GC)
  }else{
    fg_all = bind_cols(RBP_fimo_fg, 
                            TF_fimo_fg,
                            as.data.frame(fg_threemers), 
                            as.data.frame(fg_fourmers),
                            fg_GC)
    
    bg_all = bind_cols(RBP_fimo_bg, 
                       TF_fimo_bg, 
                       as.data.frame(bg_threemers), 
                       as.data.frame(bg_fourmers), 
                       bg_GC)
  }
  ### repeats
  #fg_repeats = repeatcounts(fg_regions, repeats)
  #bg_repeats = repeatcounts(bg_regions, repeats)
  
  #fg_simple_repeats = repeatcounts(fg_regions, repeats2)
  #bg_simple_repeats = repeatcounts(bg_regions, repeats2)
  
  

  return(list(fg_all, bg_all))
  
}
  
##I copied this function from https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/
vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  library(fmsb)
  
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}

write_output = function(model_list){
  
  ###write confusion matrix
  pdf(paste0("files_for_svm/consmat.pdf"), width = 7, height = 3)
  
  for(i in 1:length(model_list)){
  consmat = model_list[[i]][[3]]
  plt <- as.data.frame(consmat$table)
  plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))
  totals = plt %>% group_by(Reference) %>% summarise(tots = sum(Freq))
  toplot = merge(plt, totals)
  
  consmatout = ggplot(toplot, aes(Prediction,Reference, fill= Freq/tots)) +
    geom_tile() + geom_text(aes(label=Freq)) +
    scale_fill_gradient(low="white", high="#009194") +
    labs(x = "Reference",y = "Prediction") +
    scale_x_discrete(labels=c("foreground","background")) +
    scale_y_discrete(labels=c("background","foreground")) + 
    theme_minimal() + labs(fill='Fraction of sample')
  print(consmatout)
  }
  dev.off()
  
  ##ROC curves
  
  
  ##metrics - balanced accuracy, F1, Sensitivity, Specificity
  
  
  
  
}





