#!/usr/bin/sh

#PBS -S /usr/bin/sh
#PBS -l select=1:ncpus=4:mem=40gb
#PBS -M gauravmadappa.mandana@ieo.it
#PBS -m ea
#PBS -N merge_XXXX
#PBS -e subscripts/errXXXX.out
#PBS -o subscripts/outXXXX.out
#PBS -l walltime=2:00:00

source ~/.bashrc
source ~/.bash_profile

PATH=$PATH:/hpcnfs/data/GN2/gmandana/bin
export PATH

cd $PBS_O_WORKDIR
WD=$PBS_O_WORKDIR

samtools merge -@ 4 $WD/analysis/read2/cleanbams/nomismatch/ZZZZ_read2.nomismatch.sorted.bam \
 $WD/analysis/read2/cleanbams/nomismatch/XXXX_read2.nomismatch.sorted.bam \
 $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.sorted.bam

samtools index $WD/analysis/read2/cleanbams/nomismatch/ZZZZ_read2.nomismatch.sorted.bam

bamCoverage -p 4 --normalizeUsing RPKM --filterRNAstrand forward -bs 1 -b $WD/analysis/read2/cleanbams/nomismatch/ZZZZ_read2.nomismatch.sorted.bam -o $WD/analysis/read2bw/cleanbams/ZZZZ_read2_for.bw
bamCoverage -p 4 --normalizeUsing RPKM --filterRNAstrand reverse -bs 1 -b $WD/analysis/read2/cleanbams/nomismatch/ZZZZ_read2.nomismatch.sorted.bam -o $WD/analysis/read2bw/cleanbams/ZZZZ_read2_rev.bw

### split into fwd and rev

newdir="$WD/analysis/read2/cleanbams/nomismatch"
#
samtools view -h -Sb -@ 4 -f 16 ${newdir}/ZZZZ_read2.nomismatch.sorted.bam > ${newdir}/ZZZZ_read2.nomismatch_rev.bam
samtools view -h -Sb -@ 4 -F 16 ${newdir}/ZZZZ_read2.nomismatch.sorted.bam > ${newdir}/ZZZZ_read2.nomismatch_fwd.bam
samtools sort -@ 4 -o ${newdir}/ZZZZ_read2.nomismatch_rev.sorted.bam ${newdir}/ZZZZ_read2.nomismatch_rev.bam
samtools sort -@ 4 -o ${newdir}/ZZZZ_read2.nomismatch_fwd.sorted.bam ${newdir}/ZZZZ_read2.nomismatch_fwd.bam
samtools index ${newdir}/ZZZZ_read2.nomismatch_rev.sorted.bam
samtools index ${newdir}/ZZZZ_read2.nomismatch_fwd.sorted.bam
rm ${newdir}/ZZZZ_read2.nomismatch_rev.bam ${newdir}/ZZZZ_read2.nomismatch_fwd.bam


#### remove bams of initial sequencings
rm $WD/analysis/read2/cleanbams/nomismatch/XXXX_read2.nomismatch.sorted.bam \
 $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch.sorted.bam \
 $WD/analysis/read2/cleanbams/nomismatch/XXXX_read2.nomismatch_rev.sorted.bam \
 $WD/analysis/read2/cleanbams/nomismatch/XXXX_read2.nomismatch_fwd.sorted.bam \
 $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch_fwd.sorted.bam \
 $WD/analysis/read2/cleanbams/nomismatch/YYYY_read2.nomismatch_rev.sorted.bam 

#### remove bw of initial sequencings
rm $WD/analysis/read2bw/cleanbams/XXXX_r2_for.bw \
 $WD/analysis/read2bw/cleanbams/XXXX_r2_rev.bw \
 $WD/analysis/read2bw/cleanbams/YYYY_r2_for.bw \
 $WD/analysis/read2bw/cleanbams/YYYY_r2_rev.bw \

##### to single end

#mkdir -p ${newdir}/SE

#samtools view -H -@ 4 ${newdir}/ZZZZ_read2.nomismatch_fwd.sorted.bam > ${newdir}/ZZZZ_header_for
#samtools view -H -@ 4 ${newdir}/ZZZZ_read2.nomismatch_rev.sorted.bam > ${newdir}/ZZZZ_header_rev
#
##change the flag (fwd bam)
#samtools view -@ 4 ${newdir}/ZZZZ_read2.nomismatch_fwd.sorted.bam | awk -F'\t' -vOFS='\t' '{$2=0; print $0}' > ${newdir}/SE/ZZZZ_SE_fwd.sam
#
##change the flag (rev bam)
#samtools view -@ 4 ${newdir}/ZZZZ_read2.nomismatch_rev.sorted.bam | awk -F'\t' -vOFS='\t' '{$2=16; print $0}' > ${newdir}/SE/ZZZZ_SE_rev.sam
#
##append header
#cat ${newdir}/ZZZZ_header_for ${newdir}/SE/ZZZZ_SE_fwd.sam | samtools view -b > ${newdir}/SE/ZZZZ_SE_fwd.bam
#cat ${newdir}/ZZZZ_header_rev ${newdir}/SE/ZZZZ_SE_rev.sam | samtools view -b > ${newdir}/SE/ZZZZ_SE_rev.bam
#
#rm ${newdir}/ZZZZ_header_for ${newdir}/ZZZZ_header_rev ${newdir}/SE/ZZZZ_SE_fwd.sam ${newdir}/SE/ZZZZ_SE_rev.sam

