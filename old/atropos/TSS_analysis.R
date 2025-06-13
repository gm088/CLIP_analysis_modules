source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/ENHANCEDCLIP/feat_extr_kmer_functions.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/functions_2.R")
source("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/Shannahan/functions.R")
options(meme_bin = "/hpcnfs/home/ieo5559/meme/bin")
genome = BSgenome.Hsapiens.UCSC.hg38.masked
#this is for the coding gene annotation, remember to change the chromosome names
txdb = loadDb("/hpcnfs/data/GN2/gmandana/bin/4.1.0/home/ieo5559/insp/txdb.sqlite"); seqlevels(txdb) = paste0("chr", seqlevels(txdb))
alltx = transcripts(txdb, use.names = T, columns=c("TXNAME", "TXTYPE"))
coding = alltx[alltx$TXTYPE == "protein_coding"]

#load the Ohler
#I think I made a mistake with the strand
rev = import("/hpcnfs/data/GN2/gmandana/take2project/Ohler5primeGRO/beds/11_for.bed")
fwd = import("/hpcnfs/data/GN2/gmandana/take2project/Ohler5primeGRO/beds/11_rev.bed")

rev2 = import("/hpcnfs/data/GN2/gmandana/take2project/Ohler5primeGRO/beds/10_for.bed")
fwd2 = import("/hpcnfs/data/GN2/gmandana/take2project/Ohler5primeGRO/beds/10_rev.bed")

##merge
allreads = c(fwd, rev, fwd2, rev2)
collapsed = as.data.frame(resize(allreads, width = 1, fix = "start")) %>%
  group_by(seqnames, start, strand) %>% 
  dplyr::count("start") %>%
  dplyr::select(c("seqnames", "start", "strand", "n"))

collapsed$name = "."

TSSes = makeGRangesFromDataFrame(collapsed, 
                         start.field = "start", 
                         end.field = "start", 
                         keep.extra.columns = T)

write.table(collapsed[,c(1,2,2,5,4,3)], 
       "OhlerTSS.bed",
       quote = F,
       col.names = F,
       row.names = F,
       sep = "\t")

tesreg = subsetByOverlaps(TSSes, coding[12445])

#msa test
library(msa)
proseq = import("TSS_analysis/denovo_sorted.bed")
mysiW = customimportbed6("TSS_analysis/mysiW_fdr0.0001_gt1kb.bed")
new_siW = customimportbed6("all_siW_new.bed")
clipper_peaks = customimportbed6("TSS_analysis/clipper_new_rpkm_noucscrepeats_nosmallRNA.bed")
stratified_clipper_peaks = stratify_by_score(clipper_peaks, col_to_strat = 7)
ohlertss = customimportbed6("OhlerTSS.bed")

test = custom_homogenise(subsetByOverlaps(resize(ohlertss[ohlertss$V5 > 50], width = 1, fix = "start"), 
                                          promoters(coding, upstream = 500, downstream = 500)), ups = 500, downs = 500, "TSS")

test = custom_homogenise(subsetByOverlaps(ohlertss[ohlertss$V5 > 20], unidirectional_TU_GR), ups = 200, downs = 200, "TSS")

testseqs = getSeq(genome, resize(test, width = 100, fix = "center"))

msaobj = msa(testseqs, method = "ClustalW")
## show the results
show(msaobj)
## print the results
print(msaobj, show="alignment", showConsensus=FALSE)
## print results with custom consensus sequence
print(msaobj, show="complete", type="upperlower", thresh=c(50, 20))

#msaPrettyPrint(msaobj, output="tex", showNames="none",
#               showLogo="bottom", askForOverwrite=FALSE, verbose=FALSE,
#               showConsensus="none", alfile = "test.fasta")

MaskedAlignment <- msaobj
colM <- c(IRanges(start=1, end=20), IRanges(start=130, end=150))
colmask(MaskedAlignment) <- colM
MaskedAlignment
unmasked(MaskedAlignment)


testmat = as.matrix(testseqs)

testmat[testmat == "-"] = 0
testmat[testmat == "A"] = 1
testmat[testmat == "C"] = 2
testmat[testmat == "G"] = 2
testmat[testmat == "T"] = 1

testmat = apply(testmat, 2, as.numeric)

pheatmap(testmat, cluster_rows = F, cluster_cols = F,
         color = c("white", "green", "blue", "orange", "red"),
         breaks = c(0, 0.5, 1.5, 2.5, 3.5, 4.5))  # distances 0 to 0.5 are white, etc...

#view_motifs(consensusMatrix(msaobj))
view_motifs(consensusMatrix(testseqs))



###
writeXStringSet(getSeq(genome, promoters(resize(test, width = 1, fix = "center"), upstream = 0, downstream = 100)), 
                "test_CLIPTSS.fa")


writeXStringSet(getSeq(genome, 
                       promoters(invertStrand(resize(test, width = 1, fix = "center")), upstream = 0, downstream = 100)), 
                "test_CLIPTSS_inv_for_background.fa")









