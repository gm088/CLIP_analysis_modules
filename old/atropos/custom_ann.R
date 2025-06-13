library(rtracklayer)
library(edgeR)
library(BiocParallel)
library(GenomicRanges)
library(DESeq2)
library(dplyr)

###from my custom annotation, get the IDs of those that are affected by ZC3H4 depletion

desmat = matrix(rep(c(c(1,0,0,0),c(1,0,0,1),c(1,1,0,0),c(1,1,0,1),c(1,0,1,0),c(1,0,1,1),c(1,1,1,0),c(1,1,1,1)),3),
                nrow=4)
desmat = t(desmat)
colnames(desmat) = c("Intercept", "siCPSF3", "siCPSF3L", "AUX")

#counts = read.table("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/featcount/custom_ann/count_table.dat", header = T)

##### 10 Jul 2023 doing it for new one
counts = read.table("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/featcount/10_may_ttome/10_may_counts.dat", header = T)

dedup_counts = counts %>% 
  group_by(Geneid) %>%
  arrange(desc(Length), .by_group = T) %>%
  filter(row_number()==1) %>%
  ungroup()

filter = rowSums(dedup_counts[,c(7:30)]) > 40
counts_filt = dedup_counts[filter,]

counts_only = data.frame(counts_filt[,c(7:30)])
rownames(counts_only) = counts_filt$Geneid

y = DGEList(counts_only)
y = estimateDisp(y,desmat)
fit = glmQLFit(y, desmat, robust=T)
results = glmQLFTest(fit, coef="AUX")

column_data = t(t(colnames(counts_only)))
dds = DESeqDataSetFromMatrix(countData = counts_only, colData = column_data, design = desmat)
dds = DESeq(dds)
res = results(dds, contrast = c(0,0,0,1))

#check
all.equal(rownames(counts_only), rownames(res))
all.equal(counts_filt$Geneid, rownames(res))
table = cbind(counts_filt[,c(1,2,3,4,5)], res[,c(1,2,6)])
rownames(table) = NULL
#merge(dedup_counts, data.frame(res), by=0, all=TRUE)

####
#rpkm_scores = read.table("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/stringtie/rpkm_merged_transcripts.bed")
rpkm_scores = read.table("/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/stringtie/9May/rpkm_merged_transcripts.bed")
rownames(rpkm_scores) = rpkm_scores$V4

#different isoforms...
test = merge(rpkm_scores, table, by.x='V4', by.y='Geneid')
df = test %>% dplyr::select(-c(V1, V2, V3, V6)) %>% mutate(logpadj = -log10(padj))

#chr start stop name RPKM strand, baseMean, log2FC, logpadj
df[,c(3,4,5,1,2,6,7,8,10)]

#write.table(df[,c(3,4,5,1,2,6,7,8,10)],
#            file = "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/stringtie/rpkm_merged_transcripts_w_diffexp.bed",
#            sep = "\t",
#            quote = F,
#            row.names = F,
#            col.names = F)

write.table(df[,c(3,4,5,1,2,6,7,8,10)],
            file = "/hpcnfs/data/GN2/gmandana/take2project/t3_reseq/stringtie/9May/rpkm_merged_transcripts_w_diffexp.bed",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)





