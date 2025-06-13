library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library('biomaRt')
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/Shannahan/functions.R")

shiftStranded = function(x, value=0L,...) GenomicRanges::shift(x ,shift=value*ifelse('-'==strand(x),-1,1),...)

txdb = makeTxDbFromGFF("/hpcnfs/data/GN2/gmandana/annotation/gencode.v37.annotation.gtf")
bedfile = "/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/clipper/10Mqy_ttome/PNUTS/rep1sig_adjusted.bed"
regs = customimportbed6(bedfile)
mysiW = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/T3_ghostbusters/siW__better_annotated_20Sep.bed")
#mysiW = customimportbed6("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/DRB3/siW_2/siW_recall_6Nov.bed")


#find ones sensitive in polyA (CODING GENES)
res = readRDS("../miscellany/DEseqres.RDS")
up = rownames(res[res$padj < 0.001 & res$log2FoldChange > 0.5,])

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
#use this to find transcripts sens to zc3h4
txRes = tx2gene[tx2gene$GENEID %in% up,]$TXNAME

##

names(regs) = paste0("id", 1:length(regs))

f <- function(s) strsplit(s, ",")[[1]][1]

regs$name = NA
regs$class = NA
regs$RestrictorSuppressed = NA

maxd = 500

##start 3' UTR
threeUTRS = unlist(threeUTRsByTranscript(txdb, use.names = T))

olap1 = findOverlaps(regs, 
                     resize(threeUTRS, width = width(threeUTRS)+500, fix = "center"), 
                     ignore.strand = F)

annot1 = as.data.frame(olap1) %>% group_by(queryHits) %>% filter(row_number()==1)

regs[annot1$queryHits]$name = names(threeUTRS[annot1$subjectHits])
regs[annot1$queryHits]$class = "3.UTR"
regs[annot1$queryHits]$RestrictorSuppressed = regs[annot1$queryHits]$name %in% txRes

# more UTR
unannot = regs[is.na(regs$class)]

olap2 = findOverlaps(unannot, 
                     promoters(resize(threeUTRS, width = 1, fix = "end"), downstream = 5000, upstream = 0), 
                     ignore.strand = F)

annot2 = as.data.frame(olap2) %>% group_by(queryHits) %>% filter(row_number()==1)

unannot[annot2$queryHits]$name = names(threeUTRS[annot2$subjectHits])
unannot[annot2$queryHits]$class = "3.UTR"
unannot[annot2$queryHits]$RestrictorSuppressed = unannot[annot2$queryHits]$name %in% txRes

# promoters
unannot2 = unannot[is.na(unannot$class)]
promo = promoters(txdb, upstream = maxd*3, downstream = maxd*3)

olap3 = findOverlaps(unannot2,
                     promo, 
                     ignore.strand = F)

annot3 = as.data.frame(olap3) %>% group_by(queryHits) %>% filter(row_number()==1)
unannot2[annot3$queryHits]$name = names(promo[annot3$subjectHits])
unannot2[annot3$queryHits]$class = "promoter"
unannot2[annot3$queryHits]$RestrictorSuppressed = unannot2[annot3$queryHits]$name %in% txRes

### introns
unannot3 = unannot2[is.na(unannot2$class)]

introns = unlist(intronsByTranscript(txdb, use.names = T))

olap4 = findOverlaps(unannot3, 
                     introns, 
                     ignore.strand = F)

annot4 = as.data.frame(olap4) %>% group_by(queryHits) %>% filter(row_number()==1)

unannot3[annot4$queryHits]$name = names(introns[annot4$subjectHits])
unannot3[annot4$queryHits]$class = "intron"
unannot3[annot4$queryHits]$RestrictorSuppressed = unannot3[annot4$queryHits]$name %in% txRes

#next - paRNA
unannot4 = unannot3[is.na(unannot3$class)]
promo2 = promoters(txdb, upstream = maxd*4, downstream = maxd*2)

olap5 = findOverlaps(invertStrand(unannot4), 
                     promo2, 
                     ignore.strand = F)

annot5 = as.data.frame(olap5) %>% group_by(queryHits) %>% filter(row_number()==1)
unannot4[annot5$queryHits]$name = names(promo2[annot5$subjectHits])
unannot4[annot5$queryHits]$class = "pa-RNA"
unannot4[annot5$queryHits]$RestrictorSuppressed = unannot4[annot5$queryHits]$name %in% txRes

#??? siW
unannot5 = unannot4[is.na(unannot4$class)]

olap6 = findOverlaps(unannot5, 
                     resize(mysiW, fix = "center", width = width(mysiW)+500), 
                     ignore.strand = F)

annot6 = as.data.frame(olap6) %>% group_by(queryHits) %>% filter(row_number()==1)
unannot5[annot6$queryHits]$name = mysiW[annot6$subjectHits]$V4
unannot5[annot6$queryHits]$class = "Restrictor-sensitive"
unannot5[annot6$queryHits]$RestrictorSuppressed = T

#for now
unannot5[is.na(unannot5$name)]$name = "ambiguous"
unannot5[unannot5$name == "ambiguous"]$class = "ambiguous"
unannot5[unannot5$name == "ambiguous"]$RestrictorSuppressed = F

intermed = c(unannot5,
             unannot4[!is.na(unannot4$class)],
             unannot3[!is.na(unannot3$class)],
             unannot2[!is.na(unannot2$class)],
             unannot[!is.na(unannot$class)],
             regs[!is.na(regs$class)])

regs_annot = intermed[order(intermed$V5, decreasing = T)]
regs_annot$V9 = NULL
regs_annot$siW = countOverlaps(regs_annot, mysiW) > 0
regs_annot$RES = regs_annot$RestrictorSuppressed | regs_annot$siW

saveRDS(regs_annot, "PNUTS_rep1_peaks_annotated_23Nov.RDS")
saveRDS(as.data.frame(regs_annot[regs_annot$siW | regs_annot$class == "pa-RNA"])[,c(1,2,3,4,5,6,7,8,9,10,11)], "forDanilo.RDS")

df = as.data.frame(regs_annot)

ggplot(df, aes(class)) + 
  geom_bar(aes(fill=factor(RES)), width = 0.5, position="stack") +
  theme_minimal(base_size = 16) +
  guides(fill=guide_legend(title=""))







mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")
G_list = getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),
                values=genes,
                mart= mart)





































