if [ -z "${outdir}" ]; then echo "outdir not defined"; exit; fi

cd $outdir

#rep1
overlap_peakfi_with_bam.pl \
$rep1bam \
$inp1bam \
$rep1bed \
$rep1_readnum \
$inp1_readnum \
$outdir/rep1_norm.bed

compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl \
$outdir/rep1_norm.bed.full \
$outdir/rep1_norm_collapsed.bed \
$outdir/rep1_norm_collapsed.bed.full

#rep 2
overlap_peakfi_with_bam.pl \
$rep2bam \
$inp2bam \
$rep2bed \
$rep2_readnum \
$inp2_readnum \
$outdir/rep2_norm.bed

compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat_outputfull.pl \
$outdir/rep2_norm.bed.full \
$outdir/rep2_norm_collapsed.bed \
$outdir/rep2_norm_collapsed.bed.full

### entropy
make_informationcontent_from_peaks.pl \
$outdir/rep1_norm_collapsed.bed.full \
$rep1_readnum \
$inp1_readnum \
$outdir/rep1_norm_collapsed_entropy.bed.full \
$outdir/rep1_norm_collapsed_entropy.excessreads

make_informationcontent_from_peaks.pl \
$outdir/rep2_norm_collapsed.bed.full \
$rep2_readnum \
$inp2_readnum \
$outdir/rep2_norm_collapsed_entropy.bed.full \
$outdir/rep2_norm_collapsed_entropy.excessreads

full_to_bed.py \
--input $outdir/rep1_norm_collapsed_entropy.bed.full \
--output $outdir/rep1_norm_collapsed_entropy.bed

full_to_bed.py \
--input $outdir/rep2_norm_collapsed_entropy.bed.full \
--output $outdir/rep2_norm_collapsed_entropy.bed



