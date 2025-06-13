#!/usr/bin/sh

#PBS -S /usr/bin/sh
#PBS -l select=1:ncpus=8:mem=40gb
#PBS -M gauravmadappa.mandana@ieo.it
#PBS -m ea
#PBS -N eCLIP_YYYY
#PBS -e subscripts/errYYYY.out
#PBS -o subscripts/outYYYY.out
#PBS -l walltime=90:00:00

source ~/.bashrc
source ~/.bash_profile

PATH=$PATH:/hpcnfs/data/GN2/gmandana/bin
export PATH

cd $PBS_O_WORKDIR

ann="/hpcnfs/data/GN2/gmandana/annotation/gencode.v37.annotation.gtf"
adaptersfile="/hpcnfs/home/ieo5559/adapters_w_PCR.fa"
genomedir='/hpcnfs/data/GN2/gmandana/annotation/STAR2/genome'
genomerepeats='/hpcnfs/data/GN2/gmandana/annotation/STAR/repeats'
rawdatadir=X_RAWDATADIR_X
blacklist="/hpcnfs/data/GN2/gmandana/annotation/hg38-blacklist.v2.bed"
CLIPBL="/hpcnfs/data/GN2/gmandana/annotation/small_RNA_tothrow_and_just_tRNA_hg38_encode_BL_bed4.bed"
WD=$PBS_O_WORKDIR

##concatenate fastq files

cat $rawdatadir/XXXX/*_R1_001.fastq.gz > $WD/FASTQ/YYYY_R1.fastq.gz
cat $rawdatadir/XXXX/*_R2_001.fastq.gz > $WD/FASTQ/YYYY_R2.fastq.gz

seqkit stats $WD/FASTQ/YYYY_*.fastq.gz >> $WD/metrics/YYYY.txt
echo >> $WD/metrics/YYYY.txt

#UMI extract (read1 take the UMI+BC, read2 just UMI)

umi_tools extract -I $WD/FASTQ/YYYY_R1.fastq.gz --random-seed 1 --extract-method=regex \
--bc-pattern='^(?P<umi_1>.{12}).+' \
--bc-pattern2='^(?P<umi_2>.{10}).+' --read2-in=$WD/FASTQ/YYYY_R2.fastq.gz \
--stdout=$WD/umi/YYYY_R1_umi.fastq.gz --read2-out=$WD/umi/YYYY_R2_umi.fastq.gz

fastqc --extract -t 8 $WD/umi/YYYY_R1_umi.fastq.gz -o $WD/umi/qc
fastqc --extract -t 8 $WD/umi/YYYY_R2_umi.fastq.gz -o $WD/umi/qc

seqkit stats $WD/umi/YYYY_*.fastq.gz >> $WD/metrics/YYYY.txt
echo >> $WD/metrics/YYYY.txt

#adapter trimming

cutadapt --cores=8 --quality-cutoff 6 -m 18 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT \
-g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA \
-A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG \
-A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA \
-A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC \
-A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG \
-A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG \
-A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT \
-o $WD/trims/YYYY_R1_umi_trimmed.fastq.gz -p $WD/trims/YYYY_R2_umi_trimmed.fastq.gz \
$WD/umi/YYYY_R1_umi.fastq.gz $WD/umi/YYYY_R2_umi.fastq.gz

fastqc --extract -t 8 $WD/trims/YYYY_R1_umi_trimmed.fastq.gz -o $WD/trims/qc
fastqc --extract -t 8 $WD/trims/YYYY_R2_umi_trimmed.fastq.gz -o $WD/trims/qc

seqkit stats $WD/trims/YYYY_*_umi_trimmed.fastq.gz >> $WD/metrics/YYYY.txt
echo >> $WD/metrics/YYYY.txt

##final trimming with umi checking
conda deactivate

python3 ~/bin/excisor2.py --fname $WD/trims/YYYY_R1_umi_trimmed.fastq.gz --direction read1 \
--outname YYYY_R1_final.fastq --outdir $WD/trims --minlen 17 > metrics/YYYY_R1_excisor.out

python3 ~/bin/excisor2.py --fname $WD/trims/YYYY_R2_umi_trimmed.fastq.gz --direction read2 \
--outname YYYY_R2_final.fastq --outdir $WD/trims --minlen 17 > metrics/YYYY_R2_excisor.out

conda activate

gzip $WD/trims/YYYY_R1_final.fastq
gzip $WD/trims/YYYY_R2_final.fastq

fastqc --extract -t 8 $WD/trims/YYYY_R1_final.fastq.gz -o $WD/trims/qc
fastqc --extract -t 8 $WD/trims/YYYY_R2_final.fastq.gz -o $WD/trims/qc

seqkit stats $WD/trims/YYYY_*_final.fastq.gz >> $WD/metrics/YYYY.txt
echo >> $WD/metrics/YYYY.txt

##pair reads in case a mate survived but its mate didnt

seqkit pair -1 $WD/trims/YYYY_R1_final.fastq.gz -2 $WD/trims/YYYY_R2_final.fastq.gz -O repaired/YYYY
seqkit stats $WD/repaired/YYYY/YYYY_*_final.fastq.gz >> $WD/metrics/YYYY.txt
echo >> $WD/metrics/YYYY.txt

#map to repeat elements

mkdir $WD/repeat_mapped/YYYY

STAR \
--runThreadN 8 --genomeDir $genomerepeats \
--readFilesCommand gunzip -c \
--outFilterMultimapNmax 30 \
--outFilterMismatchNmax 10 \
--outSAMattributes All \
--outFilterMultimapScoreRange 1 \
--outSAMtype BAM Unsorted \
--outStd BAM_Unsorted \
--outReadsUnmapped Fastx \
--outFilterScoreMin 10  --outSAMattrRGline ID:foo  --alignEndsType EndToEnd \
--outFilterType BySJout \
--outFileNamePrefix $WD/repeat_mapped/YYYY/YYYY \
--readFilesIn $WD/repaired/YYYY/YYYY_R1_final.fastq.gz $WD/repaired/YYYY/YYYY_R2_final.fastq.gz > $WD/repeat_mapped/YYYY/repeat_mapped.bam 

##quality check on repeat-unmapped reads

fastqc -t 8 $WD/repeat_mapped/YYYY/YYYYUnmapped.out.mate1 -o $WD/repeat_mapped/qc
fastqc -t 8 $WD/repeat_mapped/YYYY/YYYYUnmapped.out.mate2 -o $WD/repeat_mapped/qc

##final re-pairing

mkdir $WD/ready2map
seqkit pair -1 $WD/repeat_mapped/YYYY/YYYYUnmapped.out.mate1 -2 $WD/repeat_mapped/YYYY/YYYYUnmapped.out.mate2 -O $WD/ready2map/YYYY
seqkit stats $WD/ready2map/YYYY/YYYY* >> $WD/metrics/YYYY.txt
echo >> $WD/metrics/YYYY.txt

##map to genome 

mkdir $WD/mapped/YYYY

STAR \
--runThreadN 8 --genomeDir $genomedir \
--outFilterMultimapNmax 1 \
--outSAMunmapped Within \
--outSAMattributes All \
--alignEndsType EndToEnd \
--outStd BAM_Unsorted \
--outFilterScoreMin 10 \
--outFilterType BySJout \
--outSAMtype BAM Unsorted \
--outReadsUnmapped Fastx \
--outFileNamePrefix $WD/mapped/YYYY/ \
--readFilesIn $WD/ready2map/YYYY/YYYYUnmapped.out.mate1 $WD/ready2map/YYYY/YYYYUnmapped.out.mate2 > $WD/mapped/YYYY/Aligned.bam

##generate bam and bigwig files

cd $WD/mapped/YYYY

touch redentor.txt #keep track of how much we lose at each step
cat Log.final.out | grep 'Number of input reads' >> redentor.txt

samtools flagstat Aligned.bam >> redentor.txt
echo >> redentor.txt
samtools sort -o Aligned.sorted.bam Aligned.bam
samtools index Aligned.sorted.bam

###UMI dedup
#
cd $WD/mapped/YYYY

umi_tools dedup -I Aligned.sorted.bam --log=umidedup.log --paired --extract-umi-method=read_id --method=unique > Aligned_dedup.bam

echo "umi deduplication" >> redentor.txt
samtools flagstat Aligned_dedup.bam >> redentor.txt
echo >> redentor.txt

###########remove blacklist regions
samtools sort -n -o Aligned_dedup.sorted.bam Aligned_dedup.bam

pairToBed -abam Aligned_dedup.sorted.bam -b $CLIPBL -type neither > clean_dedup.bam

echo "removed blacklisted regions" >> redentor.txt
samtools flagstat clean_dedup.bam >> redentor.txt
echo >> redentor.txt

samtools sort -o clean_dedup.sorted.bam clean_dedup.bam
samtools index clean_dedup.sorted.bam
rm clean_dedup.bam

##    take the read2 of the mapped reads, and we want to collapse to the 5' end

mkdir $WD/analysis
mkdir $WD/analysis/read2
mkdir $WD/analysis/read2bw

samtools view -hb -@ 8 -f 128 $WD/mapped/YYYY/clean_dedup.sorted.bam > $WD/analysis/read2/YYYY_read2.bam
samtools sort -@ 8 -o $WD/analysis/read2/YYYY_read2.sorted.bam $WD/analysis/read2/YYYY_read2.bam
rm $WD/analysis/read2/YYYY_read2.bam
samtools index $WD/analysis/read2/YYYY_read2.sorted.bam

## bigwig

#bamCoverage -p 8 --normalizeUsing RPKM --filterRNAstrand forward -bs 1 -b $WD/analysis/read2/YYYY_read2.sorted.bam -o $WD/analysis/read2bw/YYYY_r2_for.bw
#bamCoverage -p 8 --normalizeUsing RPKM --filterRNAstrand reverse -bs 1 -b $WD/analysis/read2/YYYY_read2.sorted.bam -o $WD/analysis/read2bw/YYYY_r2_rev.bw

## no mismatch

samtools view -H $WD/analysis/read2/YYYY_read2.sorted.bam > $WD/analysis/read2/cleanbams/nomismatch/YYYY_header
samtools view $WD/analysis/read2/YYYY_read2.sorted.bam | grep NM:i:0 > $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.sam
cat $WD/analysis/read2/cleanbams/nomismatch/YYYY_header $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.sam | samtools view -b > $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.bam

samtools sort -o $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.sorted.bam $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.bam
samtools index $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.sorted.bam
rm $WD/analysis/read2/cleanbams/nomismatch/YYYY_header $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.sam
rm $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.bam

samtools flagstat $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.sorted.bam >> $WD/mapped/YYYY/redentor.txt
echo >> $WD/mapped/YYYY/redentor.txt

bamCoverage -p 8 --normalizeUsing RPKM --filterRNAstrand forward -bs 1 -b $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.sorted.bam -o $WD/analysis/read2bw/cleanbams/YYYY_r2_for.bw
bamCoverage -p 8 --normalizeUsing RPKM --filterRNAstrand reverse -bs 1 -b $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.sorted.bam -o $WD/analysis/read2bw/cleanbams/YYYY_r2_rev.bw

### link for tracsk

ln -s $WD/analysis/read2bw/cleanbams/YYYY_r2_for.bw /hpcnfs/data/GNUH/tracks/IEO/gmandana/clip/bw/read2bw/YYYY_r2_for.bw
ln -s $WD/analysis/read2bw/cleanbams/YYYY_r2_rev.bw /hpcnfs/data/GNUH/tracks/IEO/gmandana/clip/bw/read2bw/YYYY_r2_rev.bw

### split and SE

newdir="$WD/analysis/read2/cleanbams/nomismatch"

samtools view -h -Sb -@ 4 -f 16 ${newdir}/YYYY_read2.nomismatch.sorted.bam > ${newdir}/YYYY_read2.nomismatch_rev.bam
samtools view -h -Sb -@ 4 -F 16 ${newdir}/YYYY_read2.nomismatch.sorted.bam > ${newdir}/YYYY_read2.nomismatch_fwd.bam
samtools sort -@ 4 -o ${newdir}/YYYY_read2.nomismatch_rev.sorted.bam ${newdir}/YYYY_read2.nomismatch_rev.bam
samtools sort -@ 4 -o ${newdir}/YYYY_read2.nomismatch_fwd.sorted.bam ${newdir}/YYYY_read2.nomismatch_fwd.bam
samtools index ${newdir}/YYYY_read2.nomismatch_rev.sorted.bam
samtools index ${newdir}/YYYY_read2.nomismatch_fwd.sorted.bam
rm ${newdir}/YYYY_read2.nomismatch_rev.bam ${newdir}/YYYY_read2.nomismatch_fwd.bam

##### to single end

#mkdir -p ${newdir}/SE

samtools view -H ${newdir}/YYYY_read2.nomismatch_fwd.sorted.bam > ${newdir}/YYYY_header_for
samtools view -H ${newdir}/YYYY_read2.nomismatch_rev.sorted.bam > ${newdir}/YYYY_header_rev

#change the flag (fwd bam)
samtools view ${newdir}/YYYY_read2.nomismatch_fwd.sorted.bam | awk -F'\t' -vOFS='\t' '{$2=0; print $0}' > ${newdir}/SE/YYYY_SE_fwd.sam

#change the flag (rev bam)
samtools view ${newdir}/YYYY_read2.nomismatch_rev.sorted.bam | awk -F'\t' -vOFS='\t' '{$2=16; print $0}' > ${newdir}/SE/YYYY_SE_rev.sam

#append header
cat ${newdir}/YYYY_header_for ${newdir}/SE/YYYY_SE_fwd.sam | samtools view -b > ${newdir}/SE/YYYY_SE_fwd.bam
cat ${newdir}/YYYY_header_rev ${newdir}/SE/YYYY_SE_rev.sam | samtools view -b > ${newdir}/SE/YYYY_SE_rev.bam

rm ${newdir}/YYYY_header_for ${newdir}/YYYY_header_rev ${newdir}/SE/YYYY_SE_fwd.sam ${newdir}/SE/YYYY_SE_rev.sam








