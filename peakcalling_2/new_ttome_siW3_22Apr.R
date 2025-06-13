library(dplyr)
library(stringr)
source("/Users/IEO5559/enhancedClip/Shannahan/feat_extr_kmer_functions.R")
source("~/enhancedClip/functions_2.R")
source("/Users/IEO5559/enhancedClip/Shannahan/functions.R")
### first, take bed12 of annotation and collapse and then parse
### extend 500bp upstream

annot = "/Users/IEO5559/Desktop/misc/annotations/gencode.v39.annotation.sorted.gtf"
txdb = makeTxDbFromGFF(annot)
gencod = read.table("~/Desktop/misc/annotations/gencode.v39.annotation.sorted.noheader.gtf", 
                    header = F, sep = "\t")

########### HERE make the siW gtf, with deduped names
siW_recall = "~/prog/bedfiles/mysiW/siw_3/refined/siW_22Apr_annotated.bed"

siw_ranges = customimportbed6(siW_recall)
siw_ranges_ext = resize(siw_ranges, fix = "center", width = width(siw_ranges)+2000)
siw_df = as.data.frame(siw_ranges_ext)

##dedup names ("enahncer")
deduped = vector(mode = "character", length = nrow(siw_df))
for(i in 1:nrow(siw_df)){
  
  classcounts = table(siw_df$V4)
  if(classcounts[siw_df[i,]$V4] == 1){
    #deduped[i] = siw_df[i,]$V4
    next
  }else{
    siw_df[i,]$V4 = paste0(siw_df[i,]$V4, "_", classcounts[siw_df[i,]$V4])
    #deduped[i] = paste0(siw_df[i,]$V4, "_", classcounts[siw_df[i,]$V4])
  }
  
}

siw_df$source = "CSAW"
siw_df$identifier = "transcript"
siw_df$attr = paste0("gene_id ", siw_df$V4, "; transcript_id ", siw_df$V4, "; gene_type lncRNA;")
siw_df$placeholder = "."

#siw_gtf = siw_df[,c(1,10,11,2,3,13,5,13,12)]
siw_gtf = siw_df[,c("seqnames", "source", "identifier", "start",
                    "end", "placeholder", "strand", "placeholder",
                    "attr")]
write.table(siw_gtf, file = "/Users/IEO5559/prog/DRB_TTseq/siW_ttome/siW3_22Apr_only_dedup_names.gtf",
            sep = "\t",
            quote = F, row.names = F, col.names = F)

## for combining with gencode
colnames(siw_gtf) = paste0("V", 1:dim(siw_gtf)[2])
combined = rbind(gencod, siw_gtf)

######## v2 just means the siW names are deduped
write.table(combined, file = "/Users/IEO5559/prog/DRB_TTseq/siW_ttome/4Nov_v2.gtf",
            sep = "\t",
            quote = F, row.names = F, col.names = F)

######################

## HERE MAKE THE CLIPPER FILES

######################

a = unlist(transcriptsBy(txdb, by = "gene"))
a$gene = names(a)
names(a) = NULL

reference = as.data.frame(a) %>% 
  group_by(gene) %>% 
  arrange(desc(width), .by_group = T) %>% 
  filter(row_number() == 1)

exons = unlist(exonsBy(txdb, by = "tx", use.names = T)[reference$tx_name])
exons$tx_name = names(exons)
names(exons) = NULL

premrna_length = as.data.frame(exons) %>% 
  group_by(tx_name) %>%
  dplyr::select(c("width", "tx_name")) %>% 
  summarise(mrna_length = sum(width))

####hahahahaha
final_all  = merge(reference, premrna_length, by.x = "tx_name", by.y = "tx_name") %>% 
  arrange(seqnames) %>% 
  group_by(seqnames) %>% 
  arrange((start), .by_group = T) %>% 
  ungroup()

final_exons = as.data.frame(exons) %>% arrange(seqnames) %>% 
  group_by(seqnames) %>% 
  arrange((start), .by_group = T) %>% 
  ungroup()

##use the siW above
siw_df
#siw_ranges$id = paste0("id", 1:length(siw_ranges))

siw_final = siw_df %>% 
  dplyr::select(1,2,3,4,5,6) %>%
  arrange(seqnames) %>% 
  group_by(seqnames) %>% 
  arrange((start), .by_group = T) %>% 
  ungroup()

siw_final$mrna_length = siw_final$width

colnames(siw_final) = c("seqnames", "start", "end", "width", "strand", "tx_name", "mrna_length")

combined = rbind(final_all[,c(1,2,3,4,5,6,9)], siw_final[,c(6,1,2,3,4,5,7)]) %>%
  arrange(seqnames) %>% 
  group_by(seqnames) %>% 
  arrange((start), .by_group = T) %>% 
  ungroup()

combined_exons = rbind(final_exons[,c(1,2,3,4,5,9)], siw_final[,c(1,2,3,4,5,6)]) %>% 
  arrange(seqnames) %>% 
  group_by(seqnames) %>% 
  arrange((start), .by_group = T) %>% 
  ungroup()

combined$score = 0

write.table(combined[,c(2,3,4,1,8,6)],
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F,
            file = "/Users/IEO5559/prog/DRB_TTseq/siW_ttome/4Nov.ttome_genes.bed")

combined_exons$score = 0

write.table(combined_exons[,c(1,2,3,6,7,5)],
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F,
            file = "/Users/IEO5559/prog/DRB_TTseq/siW_ttome/4Nov.ttome_exons.bed")

finalgff = combined
finalgff$source = "custom"
finalgff$type = "gene"
finalgff$score = "."
finalgff$att = "."
finalgff$summ = paste0("gene_id=", 
                       finalgff$tx_name, 
                       ";mrna_length=", 
                       finalgff$mrna_length, 
                       ";premrna_length=", 
                       finalgff$width)

write.table(finalgff[,c("seqnames", "source", "type", "start", "end", "score", "strand", "att", "summ")],
            quote = F,
            sep = "\t", 
            row.names = F,
            col.names = F,
            file = "/Users/IEO5559/prog/DRB_TTseq/siW_ttome/4Nov.ttome.AS.STRUCTURE.COMPILED.gff")








################ define first 5kb (for sense-antisense, featurecounts)

siw_df_5kb = as.data.frame(promoters(siw_ranges, upstream = 0, downstream = 2000))
siw_df_5kb$source = "CSAW"
siw_df_5kb$identifier = "transcript"
siw_df_5kb$attr = paste0("gene_id ", siw_df_5kb$V4, "; gene_type lncRNA;")
siw_df_5kb$placeholder = "."

siw_gtf_5kb = siw_df_5kb[,c(1,10,11,2,3,13,5,13,12)]
colnames(siw_gtf_5kb) = paste0("V", 1:dim(siw_gtf_5kb)[2])

write.table(siw_gtf_5kb, file = "/Users/IEO5559/prog/DRB_TTseq/7Nov24_ttome/siW_2kb.gtf",
            sep = "\t",
            quote = F, row.names = F, col.names = F)
