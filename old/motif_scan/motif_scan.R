library(dplyr)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(fitdistrplus)
library(GenomicRanges)
library(GenomicAlignments)
library(slider)
library(csaw)
library(ggplot2)
library(stringr)
library(universalmotif)
source("functions_2.R")
library(kmer)
library(seqLogo)
library(ggseqlogo)
library(reshape2)
library(pheatmap)
library(dendextend)

genome = BSgenome.Hsapiens.UCSC.hg38.masked

fimo = read.table("~/prog/clip/IP7_and_DCIP/13_Apr_CSTF/all_motifs_tartaglia_5kbwin_centred/fimo.tsv", 
                  header = T, sep = "\t")
fimo$protein = unlist(lapply(str_split(fimo$motif_id, "\\|"), function(x){x[2]}))

motif_file = "/Users/IEO5559/Downloads/Homo_sapiens/1_more.meme"
motifs = read_meme(motif_file)
names = lapply(motifs, function(x){str_split(x["name"], "\\|")[[1]][2]})
ids = lapply(motifs, function(x){x["name"]})

#can't amke matrix cos rownames are not identical
for(i in 1:length(motifs)){
  motifs[[i]]["name"] = names[[i]]
}

##how many hits per protein?
hits_per_protein = as.data.frame(fimo %>% group_by(protein) %>% tally())
rownames(hits_per_protein) = hits_per_protein$protein

###pick the motif that represents the most number of hits
##for each protein, is it dominated by one motif?
fimo %>% group_by(protein) %>% group_by(motif_id, .add = T) %>% tally()

fimo %>% 
  group_by(protein) %>% 
  group_by(motif_id, .add = T) %>% 
  tally() %>% 
  mutate(uniqid = gsub("\\|", "", motif_id)) %>%
  ungroup() %>% ungroup() %>%
  ggplot(., aes(x = as.factor(protein), y = n, fill = uniqid)) + 
  geom_bar(stat = "identity", position="stack") + 
  theme(legend.position="none")

top_motifs_by_protein = fimo %>% 
  group_by(protein) %>% 
  group_by(motif_id, .add = T) %>% 
  tally() %>% 
  arrange(desc(n), .by_group = T) %>%
  filter(row_number() == 1)

rownames(top_motifs_by_protein) = top_motifs_by_protein$protein

####how many peaks per motif?
peaks_per_motif = fimo[fimo$motif_id %in% top_motifs_by_protein$motif_id,] %>% group_by(protein) %>% 
  summarise(uniq = n_distinct(sequence_name))
rownames(peaks_per_motif) = peaks_per_motif$protein

#I'm making a matrix of top motoifs per prot

top_motifs = motifs[unlist(lapply(motifs, function(x){x["name"] %in% top_motifs_by_protein$motif_id}))]
motmat = as.matrix(top_motifs)
rownames(motmat) = unlist(lapply(top_motifs, function(x){str_split(x["name"], "\\|")[[1]][2]}))

#can use name, consensus, motif(forPWM), etc
motifs[[1]]["name"]
view_motifs(motifs[10:20])

ggseqlogo(lapply(1:10, function(x){motifs[[x]]["motif"]}), ncol = 1)

view_motifs(motifs[unlist(lapply(motifs, function(x){x["name"] == "CBX7"}))])

##single line heatmap
par(mfrow = c(2,1))
matplot(colMeans(matrix), type = "l")

check = hist(fimo_3$start, breaks=500)
pheatmap(t(as.matrix(check$count)), cluster_rows = F, cluster_cols = F,
         labels_col = -2.5:2.5, 
         color = c(colorRampPalette(c("purple", "yellow", "orange"))(50)), 
         breaks = c(seq(0,10, length.out = 51)))


plot(density(fimo[fimo$protein == "CBX7",]$start))

####make a matrix, with proteins as rows, and the density/histogram of start 

hists = fimo[fimo$motif_id %in% top_motifs_by_protein$motif_id,] %>% 
  dplyr::select(protein, start) %>% 
  tidyr::pivot_wider(names_from = protein, values_from = start)

dens = lapply(t(hists), density)

densmat = t(sapply(dens, function(x){x[[2]]}))
rownames(densmat) = colnames(hists)

pheatmap(densmat, cluster_rows = T, cluster_cols = F,
         labels_col = -2.5:2.5)

##check this
pheatmap(densmat[140:160,], cluster_rows = F, cluster_cols = F,
         labels_col = -2.5:2.5)

view_motifs(motifs[unlist(lapply(motifs, function(x){x["name"] == "RBM4B"}))],
            tryRC = F)

###

test = kmeans(densmat, centers = 4)

pheatmap(densmat[test$cluster == 3,], cluster_rows = F, cluster_cols = F,
         labels_col = -2.5:2.5,
         color = c(colorRampPalette(c("white", "black"))(50)), 
         breaks = c(seq(0,0.0008, length.out = 51)))

ggseqlogo(lapply(motmat[names(test$cluster)[(test$cluster == 1)],], function(x){x["motif"]}), ncol = 4) + 
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        text = element_text(size = 5))

plotclust = function(k){
  
  hpp_subset = top_motifs_by_protein[names(test$cluster)[(test$cluster == k)],]
  ppm_subset = peaks_per_motif[names(test$cluster)[(test$cluster == k)],]
  rowlabels = paste0("N=",
                     ppm_subset$uniq,
                     "   ",
                     hpp_subset$protein, 
                     ", n=", 
                     hpp_subset$n)
  pdf(paste("motif_scan/", "dist_clust", k, ".pdf", sep = ""), width = 14, height = 14)
  pheatmap(densmat[test$cluster == k,], cluster_rows = F, cluster_cols = F,
           labels_col = -2.5:2.5,
           labels_row = rowlabels)
  dev.off()
  
  #for pdf
  factor = length(names(test$cluster)[(test$cluster == k)])
  
  pdf(paste("motif_scan/", "motif_clust", k, ".pdf", sep = ""), width = 10, height = 2*factor)
  subset = motmat[names(test$cluster)[(test$cluster == k)],]
  names = unlist(lapply(subset, function(x){x["name"]}))
  mots = lapply(subset, function(x){x["motif"]})
  print(ggseqlogo(mots, ncol = 1) + 
    theme(axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          text = element_text(size = 10)))
  #print(motif_plot)
  dev.off()
}

plotclust(4)


#####sources of motif
unlist(lapply(
  str_split(unlist(
  lapply(str_split(top_motifs_by_protein$motif_id, "\\|"), function(x){x[6]})
  ), ","), function(x){x[1]}
))

thing = table(unlist(lapply(
  str_split(unlist(
    lapply(str_split(top_motifs_by_protein$motif_id, "\\|"), function(x){x[6]})
  ), ","), function(x){x[1]}
))) %>% as.data.frame() %>% arrange(desc(Freq)) %>% filter(Freq > 1)

x = barplot(thing$Freq, xaxt="n")
labs <- thing$Var1
text(cex=0.9, x=x-0.8, y=-8, labs, xpd=TRUE, srt=45)






