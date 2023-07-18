source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/Shannahan/functions.R")
options(meme_bin = "/hpcnfs/home/ieo5559/meme/bin")
genome = BSgenome.Hsapiens.UCSC.hg38.masked
#this is for the coding gene annotation, remember to change the chromosome names
txdb = loadDb("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/txdb.sqlite")
seqlevels(txdb) = paste0("chr", seqlevels(txdb))
alltx = transcripts(txdb, use.names = T, columns=columns(txdb))
coding = alltx[alltx$TXTYPE == "protein_coding"]
proseq = import("TSS_analysis/denovo_sorted.bed")
mysiW = customimportbed6("TSS_analysis/mysiW_fdr0.0001_gt1kb.bed")
clipper_peaks = customimportbed6("TSS_analysis/clipper_new_rpkm_noucscrepeats_nosmallRNA.bed")
new_siw = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/new_siW/all_siW_new.bed")
RNA_ttome = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/Shannahan/ALLTSS_newttome_1.bed")
Pol_ttome = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/Shannahan/denovo_promprox_new.bed")
colnames(mcols(RNA_ttome)) = c("txname", "rpkm", "basemeans", "log2FC", "pval")
colnames(mcols(Pol_ttome)) = c("txname", "Pol_density")
NDR = reduce(customimportbed6("/hpcnfs/data/GN2/gmandana/annotation/ATAC_HCT116_peaks.bed"), 
             min.gapwidth = 500)
NDR$id = paste0("NDR_", 1:length(NDR))
ohlertss = customimportbed6("OhlerTSS.bed")

### first, generate all system TUs; non redundant
augmented = ttome_augment(RNA_ttome = RNA_ttome, 
                               Pol_ttome = Pol_ttome, 
                               NDR = NDR, 
                               ext_for_olap = 500)

augmented$NDR = sapply(str_split(names(augmented), "\\."), function(x){x[1]})
names(augmented) = NULL
#write.table(as.data.frame(augmented)[,c(1,2,3,6,7,5,8,9,10,11,12)],
#            sep = "\t",
#            quote = F,
#            row.names = F,
#            col.names = F,
#            file = "augmented_NDR500bpmerged.bed")

systemTUs = augmented
systemTUs = readRDS("augmented_15Jul.RDS")

#tagging tx per NDR for pairing later
systemTUs = makeGRangesFromDataFrame(as.data.frame(systemTUs) %>% 
                                       group_by(NDR) %>% 
                                       mutate(NDR_tag = paste0(NDR, "_tx", row_number())),
                                     keep.extra.columns = T)

####recount here

#### recount_augmented.R
#saveRDS(augmented_new_ranges, "augmented_16Jul_recounted.RDS")

#systemTUs = augmented_new_ranges
####

systemTUs = readRDS("augmented_16Jul_recounted.RDS")

##### now construct 
# first, overlap with pcg to get the protein coding TUs
systemTUs$iscoding = countOverlaps(systemTUs, coding) > 0

##confidence at the cost of sensitivity
systemTUs = systemTUs[!is.na(systemTUs$logCPM)]
systemTUs = systemTUs[is.na(systemTUs$alt)]

###get Ohler ####

systemTUs_w_ohler = unlist(GRangesList(
  lapply(1:length(systemTUs), function(x){get_5p_GRO(ohlertss, NDR, systemTUs[x])})
  ))


####first, bidirectional (pcg-pcg)
# same NDR, both coding
bidir = makeGRangesFromDataFrame(as.data.frame(systemTUs) %>% 
                           dplyr::filter(iscoding == T) %>% 
                           group_by(NDR) %>%  ###operate by NDR
                           dplyr::filter(n() > 1) %>%   ###must be > 1 tx per NDR
                           dplyr::filter(length(unique(strand)) > 1),   ##must be diff strand
                         keep.extra.columns = T)

####second, divergent
# same NDR, must have both coding and noncoding
div = makeGRangesFromDataFrame(as.data.frame(systemTUs) %>% 
                                   group_by(NDR) %>%  ###operate by NDR
                                   dplyr::filter(n() > 1) %>%   ###must be > 1 tx per NDR
                                   dplyr::filter(length(unique(strand)) > 1) %>%   ##must be diff strand
                                   dplyr::filter(length(unique(iscoding)) > 1),   ##one must be noncoding
                                 keep.extra.columns = T)


####third, unidirectional - don't trust this, it suffers from thresholding
# only one tx for a given NDR
unidir = makeGRangesFromDataFrame(as.data.frame(systemTUs) %>% 
                                 group_by(NDR) %>%  ###operate by NDR
                                 dplyr::filter(n() == 1),   ###must be 1 tx per NDR
                               keep.extra.columns = T)

##enhancer
enh = makeGRangesFromDataFrame(as.data.frame(systemTUs) %>%
                                 dplyr::filter(is.na(alt)) %>% 
                                 dplyr::filter(iscoding == F) %>% 
                                 group_by(NDR) %>%
                                 dplyr::filter(n() > 1) %>%
                                 dplyr::filter(length(unique(strand)) > 1),
                               keep.extra.columns = T)

length(bidir) + length(unidir) + length(div) + length(enh)

length(systemTUs)

table(sapply(div$ohler, length) == 0)


div_sensitive = subsetByOverlaps(div, new_siw)

div_sensitive_mates = check = div[div$NDR %in% div_sensitive$NDR & !(div$NDR_tag %in% div_sensitive$NDR_tag)]

as.data.frame(div) %>%
  ggplot(aes(x = iscoding, y = logFC)) + 
  geom_boxplot(aes(color=factor(iscoding)), outlier.shape=NA, width = 0.5, lwd=1.25, notch = T) + 
  theme_minimal()

#mate-wise
check = as.data.frame(div) %>% 
  group_by(NDR) %>% 
  arrange(iscoding, .by_group = T) %>% 
  mutate(diff = abs(lag(logFC)/logFC))




