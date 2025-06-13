if [ -z "${outdir}" ]; then echo "outdir not defined"; exit; fi

cd $outdir

parse_idr_peaks.pl \
$outdir/rep1vrep2.idr.out \
$outdir/rep1_norm_collapsed_entropy.bed.full \
$outdir/rep2_norm_collapsed_entropy.bed.full \
$outdir/rep1vrep2.idr.bed

#second norm
overlap_peakfi_with_bam.pl \
$rep1bam \
$inp1bam \
$outdir/rep1vrep2.idr.bed \
$rep1_readnum \
$inp1_readnum \
$outdir/rep1vrep2.idr_merged_1.bed

overlap_peakfi_with_bam.pl \
$rep2bam \
$inp2bam \
$outdir/rep1vrep2.idr.bed \
$rep2_readnum \
$inp2_readnum \
$outdir/rep1vrep2.idr_merged_2.bed

#final
get_reproducing_peaks.pl \
$outdir/rep1vrep2.idr_merged_1.bed.full \
$outdir/rep1vrep2.idr_merged_2.bed.full \
$outdir/reproducible_peaks.01.bed.full \
$outdir/reproducible_peaks.02.bed.full \
$outdir/reproducible_peaks.bed \
$outdir/reproducible_peaks.custombed \
$outdir/rep1_norm_collapsed_entropy.bed.full \
$outdir/rep2_norm_collapsed_entropy.bed.full \
$outdir/rep1vrep2.idr.out

##sort

sort -k1,1 -k2,2n $outdir/reproducible_peaks.bed > $outdir/reproducible_peaks_sorted.bed





