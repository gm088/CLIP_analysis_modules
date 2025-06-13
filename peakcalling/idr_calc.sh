if [ -z "${outdir}" ]; then echo "outdir not defined"; exit; fi

source ~/.bashrc
source ~/.bash_profile

cd $outdir

idr \
--samples \
$outdir/rep1_norm_collapsed_entropy.bed \
$outdir/rep2_norm_collapsed_entropy.bed \
--input-file-type bed \
--rank 5 \
--peak-merge-method max \
--output-file $outdir/rep1vrep2.idr.out



