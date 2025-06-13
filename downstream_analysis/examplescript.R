
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("IRanges")

## rtracklayer allows importing and exporting of many useful filetypes
## such as bed files, bigwig files, etc.
BiocManager::install("rtracklayer")

## r tracklayer automatically stores your bed corrdinates as a granges object
mybedranges = rtracklayer::import("path/to/your/bedfile")

## you can also construct a Granges object manually...
range1 = GRanges(seqnames = "chr2", 
                 ranges = IRanges(start = 231676977,
                                  end = 231695127),
                 strand = "+")

range2 = GRanges(seqnames = "chr2", 
                 ranges = IRanges(start = 231666323,
                                  end = 231684473),
                 strand = "+")

####### these functions are very useful - take a look at the descriptions

?`findOverlaps,GenomicRanges,GenomicRanges-method`

#######

##returns range1
subsetByOverlaps(range1, range2, ignore.strand = F)

##returns range2
subsetByOverlaps(range2, range1, ignore.strand = F)

##returns nothing
subsetByOverlaps(range2, range1, minoverlap = (width(range1)*0.5))

##concatenate and merge
bothranges = c(range1, range2)
merged = reduce(bothranges)





