source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/Shannahan/functions.R")
library(hues)

#load the annotation
mane = makeTxDbFromGFF("/hpcnfs/data/GN2/gmandana/annotation/MANE.GRCh38.v1.0.ensembl_genomic.gtf")
pctx = transcripts(mane, use.names = T)

procap = import("~/unifiedCLIP/PROCAP_maxima_rep1_v2.bed")

#the peaksfile
bedfile = "/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_26Oct/reproducible_peaks_reseq_rpkm.bed"
regs = customimportbed6(bedfile)
regs[1]$V5 = 1; regs[1]$V7 = 0.1; regs[1]$V8 = 3; regs[1]$V9 = 0

#mysiW = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/T3_ghostbusters/nctx_refinement/siW_refined_annotated_fantom.bed")
mysiW = readRDS("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/T3_ghostbusters/ZC3H4_supp_tx_PROCAP_FANTOM.RDS")

#clusters = makeGRangesFromDataFrame(
#  readRDS("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/peak_files_26Oct/clusters_final.RDS"),
#  keep.extra.columns = T)
clusters = makeGRangesFromDataFrame(
  readRDS("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/ZC3H4_final/clusters_final.RDS"),
  keep.extra.columns = T)

##### for the sake of the motif search, use regs - expand 100nt either side FOR THE MOTIF SEARCH,
##### append all info

names(regs) = 1:length(regs)
regsmerged = csaw::mergeWindows(resize(regs, width = width(regs) + 200, fix = "center"), 
                                tol = 0, ignore.strand = F)
regs$clus = regsmerged$ids
regsmerged$regions$clus = regsmerged$ids

sumtab = as.data.frame(regs) %>% 
  group_by(clus) %>% 
  dplyr::select(-c("seqnames", "strand", "start", "end", "width", "V9")) %>% 
  summarise_all(paste, collapse = ",")

meanscores = as.data.frame(regs) %>% 
  group_by(clus) %>% 
  summarise(mean_score = mean(V5))

toappend = as_tibble(merge(sumtab, meanscores))

opregs = regsmerged$regions
mcols(opregs) = toappend
names(opregs) = 1:length(opregs)

###### motif information
mots = switch_alph(read_meme("../miscellany/streme_final.txt"))
motifs = c(mots[[1]], mots[[5]], mots[[12]])
fimo_bg_generator(opregs, "fimo_bfile_opregs")

scan = FIMO_parse_coords_stranded(regions = opregs, thresh = 0.001,
                                  motif = motifs, fimobfile = "fimo_bfile_opregs",
                                  countseq = T)

opreghits = runFimo(getSeq(BSgenome.Hsapiens.UCSC.hg38.masked, opregs), 
                    thresh = 0.001,
                    motif = motifs, 
                    parse_genomic_coord = F,
                    bfile = "fimo_bfile_opregs", 
                    norc = T, meme_path = "/hpcnfs/home/ieo5559/meme/bin/")

#### collect and append motif information
nhits = as.data.frame(table(as.data.frame(opreghits)$seqnames))
motidata = as.data.frame(opreghits) %>% group_by(seqnames) %>% 
  dplyr::select(c(seqnames, motif_id, score)) %>% 
  summarise_all(paste, collapse = ";")

summotidata = as.data.frame(opreghits) %>% group_by(seqnames) %>% 
  dplyr::select(c(seqnames, score)) %>% 
  mutate(summotscore = sum(score)) %>%
  dplyr::select(c(seqnames, summotscore))

toappend2 = as_tibble(merge(merge(motidata, nhits, by.x = "seqnames", by.y = "Var1"),
                            summotidata), by.x = "seqnames", by.y = "seqnames")


#initialise metadata columns
opregs$nhits = NA
opregs$motscores = NA
opregs$summotscores = NA

opregs[toappend2$seqnames]$nhits = toappend2$Freq
opregs[toappend2$seqnames]$motscores = toappend2$score
opregs[toappend2$seqnames]$summotscores = toappend2$summotscore

######### final dataframe

opregs = opregs[order(opregs$mean_score, decreasing = T)]

########

############## sample A1 - do peaks with more hits have a higher score?

opregs_mut1 = as.data.frame(opregs) %>% 
  dplyr::select(c(mean_score, nhits, summotscores)) %>% 
  mutate(nhits2 = case_when(is.na(nhits) ~ 0,
                            nhits == 1 ~ 1,
                            nhits == 2 ~ 2,
                            nhits == 3 ~ 3,
                            nhits > 3 ~ 4)) %>%
  mutate_at("nhits2", as.factor)

ggplot(opregs_mut1, aes(nhits2, mean_score, fill = nhits2)) + 
  geom_violin(width=0.6, alpha = 0.4, linewidth = 1.0) + 
  geom_boxplot(width=0.1, color="grey", alpha=0.8, 
               position=position_dodge(0.6), linewidth = 1.0) + 
  scale_x_discrete(breaks=c(0.45,1,2,3,4),
                   labels=paste0("n=", table(opregs_mut1$nhits2))) +
  ylim(c(0, 60)) + 
  scale_fill_manual(values=cbbPalette) + 
  theme_linedraw()



##### first tag if its PROCAP-proximal
subsetByOverlaps(promoters(opregs, upstream = 2000, downstream = 500), procap)
procap_assoc = (as.data.frame(findOverlaps(promoters(opregs, upstream = 2000, downstream = 500), procap)) %>% 
  group_by(queryHits) %>% 
  filter(row_number() == 1))$queryHits

opregs$procap_assoc = F
opregs[procap_assoc]$procap_assoc = T

### split by pcg/not pcg

opregs$coding = F
pcg_inds = (as.data.frame(findOverlaps(opregs, pctx)) %>% 
                  group_by(queryHits) %>% 
                  filter(row_number() == 1))$queryHits

opregs[pcg_inds]$coding = T


#### bar plot splits

############## sample B1 - do score of coding/noncoding/TSS/nonTSS
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
swatch(cbbPalette)

############## FIG 2D

as.data.frame(opregs) %>%
  dplyr::select(c(procap_assoc, nhits, coding)) %>% 
  mutate(hashit = !is.na(nhits)) %>% 
  dplyr::select(c(hashit, coding, procap_assoc)) %>% 
  group_by(hashit, coding, procap_assoc) %>% 
  tally() %>%
  ggplot(aes(procap_assoc, n, fill = hashit)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~coding) + 
  theme_classic() + 
  scale_fill_manual(values=cbbPalette)

as.data.frame(opregs) %>%
  dplyr::select(c(procap_assoc, coding, siW)) %>% 
  mutate(overlap = !is.na(siW)) %>% 
  group_by(procap_assoc, coding, overlap) %>%
  tally() %>%
  ggplot(aes(procap_assoc, n, fill = overlap)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~coding) + 
  theme_classic() + 
  scale_fill_manual(values=cbbPalette)

################


#counts = RNAcount(promoters(opregs, upstream = 0, downstream = 1000), factor = 1)
#save.image("eCLIP_characterisation_jul24.RData")

############ add the actual siW for the final excel




opregs$siW = NA
olapobj = findOverlaps(opregs, mysiW)

siW_inds = (as.data.frame(olapobj) %>% 
              group_by(queryHits) %>% 
              filter(row_number() == 1))$queryHits

siWolap = mysiW[subjectHits(olapobj)]
toapp = paste0(as.vector(seqnames(siWolap)), ":", start(siWolap), "-", end(siWolap), ":", as.vector(strand(siWolap)))

#opregs[queryHits(olapobj)]$siW = mysiW[subjectHits(olapobj)]$V4
opregs[queryHits(olapobj)]$siW = toapp

saveRDS(opregs, "~/ENHANCEDCLIP/ZC3H4_final/v2/opregs.RDS")





##########################

#### coming back for RNAseq counts

opregs = readRDS("characterisation_16Jul24/opregs.RDS")
opregs$issiW = !is.na(opregs$siW)

opregs$controlcounts = counts[,1]
opregs$depZcounts = counts[,2]

as.data.frame(opregs) %>%
  dplyr::select(c(mean_score, coding, issiW, controlcounts, depZcounts)) %>% 
  mutate_at("coding", as.factor) %>%
  mutate_at("issiW", as.factor) %>%
  pivot_longer(., cols = ends_with("counts"), names_repair = "unique") %>%
  mutate_at("name", as.factor) %>%
  mutate(name = fct_relevel(name, "controlcounts", "depZcounts")) %>%
  ggplot(aes(issiW, log2(value*1e9))) + 
  geom_violin(aes(fill = name), width=0.6, alpha = 0.4, linewidth = 1.0) + 
  geom_boxplot(aes(fill = name), width=0.1, color="grey", alpha=0.8, 
               position=position_dodge(0.6), linewidth = 1.0) + 
  theme_classic() + 
  scale_fill_manual(values=cbbPalette) + 
  facet_wrap(~coding)

### p-values

yyy = opregs[opregs$coding == F & opregs$issiW == T]

t.test(yyy$depZcounts,
            yyy$controlcounts)$p.value

wilcox.test(log2(yyy$depZcounts),
            log2(yyy$controlcounts))$p.value











