source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/Shannahan/functions.R")

#load the annotation
txdb = loadDb("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/txdb.sqlite")
seqlevels(txdb) = paste0("chr", seqlevels(txdb))
alltx = transcripts(txdb, use.names = T, columns = c('TXTYPE', 'CDSNAME'))
threeutr = unlist(threeUTRsByTranscript(txdb))
procap = import("PROCAP_maxima_rep1_v2.bed")

#the peaksfile
bedfile = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_26Oct/reproducible_peaks_reseq_rpkm.bed"
regs = customimportbed6(bedfile)
mysiW = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/T3_ghostbusters/siW__better_annotated_20Sep.bed")
clusters = makeGRangesFromDataFrame(
  readRDS("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_26Oct/clusters_final.RDS"),
  keep.extra.columns = T)

#### get pcg
pcg = alltx[alltx$TXTYPE == 'protein_coding']

##### first tag if its PROCAP-proximal
subsetByOverlaps(promoters(clusters, upstream = 2000, downstream = 500), procap)
procap_assoc = (as.data.frame(findOverlaps(promoters(clusters, upstream = 2000, downstream = 500), procap)) %>% 
  group_by(queryHits) %>% 
  filter(row_number() == 1))$queryHits

clusters$procap_assoc = F
clusters[procap_assoc]$procap_assoc = T

### split by pcg/not pcg

clusters$coding = F
pcg_inds = (as.data.frame(findOverlaps(clusters, pcg)) %>% 
                  group_by(queryHits) %>% 
                  filter(row_number() == 1))$queryHits

clusters[pcg_inds]$coding = T

#########coding subcategorising
clus_pcg = subsetByOverlaps(clusters, pcg)

#promoters
clus_prom = subsetByOverlaps(clus_pcg, promoters(pcg, upstream = 1000, downstream = 2000))
p1 = subsetByOverlaps(clus_pcg, promoters(pcg, upstream = 1000, downstream = 2000), invert = T)


clus_utr = subsetByOverlaps(clus_pcg, threeutr)
p2 = subsetByOverlaps(clus_pcg, threeutr, invert = T)


####noncoding
#how many unnannotated here?
clus_notpcg = subsetByOverlaps(clusters, pcg, invert = T)
#clus_notpcg[is.na(clus_notpcg$overlap)]
subsetByOverlaps(alltx, clus_notpcg[is.na(clus_notpcg$overlap)])



#### possible visualisations

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

as.data.frame(clusters)[,c(12:19)] %>% 
  pivot_longer(., cols = ends_with("counts"), names_repair = "unique") %>%
  mutate_at("name", as.factor) %>%
  mutate(name = fct_relevel(name, "ctlcounts", "auxcounts")) %>%
  mutate("siW" = as.factor(!is.na(overlap))) %>%
  ggplot(aes(coding, log2(value*1e6))) + 
  geom_violin(aes(fill = name), width=0.6, alpha = 0.4, linewidth = 1.0) + 
  geom_boxplot(aes(fill = name), width=0.1, color="grey", alpha=0.8, 
               position=position_dodge(0.6), linewidth = 1.0) + 
  theme_classic() + 
  scale_fill_manual(values=cbbPalette) + 
  xlab("")

as.data.frame(clusters)[,c(12,18,19)] %>%
  mutate("siW" = as.factor(!is.na(overlap))) %>%
  group_by(procap_assoc, coding, siW) %>%
  summarise(n = n()) %>%
  ggplot(., aes(x = procap_assoc, y = n, fill=siW)) + 
  geom_bar(stat = "identity", width=0.5, position = "stack") + 
  facet_wrap(~ coding, ncol = 2) +
  theme_classic() + scale_fill_manual(values=cbbPalette)

wilcox.test(log2(clusters[clusters$siW == T,]$ctlcounts+1), 
              log2(clusters[clusters$siW == T,]$auxcounts+1))$p.val

wilcox.test(log2(clusters[clusters$siW == F,]$ctlcounts+1), 
            log2(clusters[clusters$siW == F,]$auxcounts+1))$p.val



























