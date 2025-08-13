# Custom analysis pipeline for enhanced Cross-linking and Immunoprecipitation (eCLIP) analysis

Running enhanced CLIP analysis on HPC cluster with PBS Pro

## Description

This repo contains the scripts needed from start to finish

* preprocessing - UMI extraction, adapter trimming, mapping to human genome
* peakcalling - generate custom annotation for peakcaller, call peaks using CLIPpers [IDR](https://www.encodeproject.org/pipelines/ENCPL357ADL/)
* peakcalling_2 - calibrate scores of called peaks for final ranking, de-novo motif search

The peak calling models reads on a transcripts to a poisson distribution in order to call peaks, with parameter equal to read count per unit length - see here. 

The score calibration adjusts the score of a peak (the -log10 p-value) by the abundance of the underlying transcript. The former is simply divided by the latter. Though perhaps overly simple, the results are the same if the transcript abundances in the system are fitted to a normal distribution, and the fitted values used for the calibration.

```
example
```

## Misc

* [Brian Yee's CLIPper docker images](https://hub.docker.com/r/brianyee/clipper/) - Used for peak calling
* [Brian Yee's IDR docker images](https://hub.docker.com/r/brianyee/merge_peaks/) - Used for calculation of enrichment and IDR
