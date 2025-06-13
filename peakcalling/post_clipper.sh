#!/usr/bin/sh

#PBS -S /usr/bin/sh
#PBS -l select=1:ncpus=8:mem=30gb
#PBS -M gauravmadappa.mandana@ieo.it
#PBS -m ea
#PBS -N post_clipper
#PBS -e err.out
#PBS -o out.out
#PBS -l walltime=20:00:00

source ~/.bashrc
source ~/.bash_profile
PATH=$PATH:/hpcnfs/data/GN2/gmandana/bin
export PATH

cd $PBS_O_WORKDIR
WD=$PBS_O_WORKDIR

#intialise variables here
#these are substituted by variables defined in the generator script
outdir=${WD}/OUTDIR
rep1bam=REP1BAM
rep1bed=REP1BED
rep2bam=REP2BAM
rep2bed=REP2BED

#num. reads
samtools view -c ${rep1bam} > rep1_readnum.txt
samtools view -c ${rep2bam} > rep2_readnum.txt
rep1_readnum="${WD}/rep1_readnum.txt"
rep2_readnum="${WD}/rep2_readnum.txt"

#leave this hardcoded for now
inp1bam="/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/13_read2.nomismatch.sorted.bam"
inp2bam="/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/read2/cleanbams/nomismatch/13_read2.nomismatch.sorted.bam"
inp1bed="/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/clipper/10Mqy_ttome/ZC3H4_old/13_read2_nomismatch.bed"
inp2bed="/hpcnfs/data/GN2/gmandana/take2project/clip_2_electric_boogaloo/analysis/clipper/10Mqy_ttome/ZC3H4_old/13_read2_nomismatch.bed"

#inp num. reads
samtools view -c ${inp1bam} > inp1_readnum.txt
samtools view -c ${inp2bam} > inp2_readnum.txt
inp1_readnum="${WD}/inp1_readnum.txt"
inp2_readnum="${WD}/inp2_readnum.txt"

###singularity images
prefix="/hpcnfs/data/GN2/gmandana/singularity_images"
mergepeaks="${prefix}/merge_peaks_0.1.0.sif"
idr="${prefix}/idr_2.0.2.sif"
python="${prefix}/eclip_0.7.0_python.sif"

###scripts
scriptsdir="/hpcnfs/data/GN2/gmandana/bin/clipper_idr_scripts"
norm_compress_entropy_script="${scriptsdir}/first_norm.sh"
idr_script="${scriptsdir}/idr_calc.sh"
norm_using_idr_peaks_script="${scriptsdir}/second_norm.sh"

mkdir -p $outdir

SINGULARITYENV_outdir=${outdir} SINGULARITYENV_rep1bam=${rep1bam} SINGULARITYENV_rep1bed=${rep1bed} \
SINGULARITYENV_rep2bam=${rep2bam} SINGULARITYENV_rep2bed=${rep2bed} SINGULARITYENV_inp1bam=${inp1bam} \
SINGULARITYENV_inp2bam=${inp2bam} SINGULARITYENV_inp1bed=${inp1bed} SINGULARITYENV_inp2bed=${inp2bed} \
SINGULARITYENV_inp1_readnum=${inp1_readnum} SINGULARITYENV_inp2_readnum=${inp2_readnum} \
SINGULARITYENV_rep1_readnum=${rep1_readnum} SINGULARITYENV_rep2_readnum=${rep2_readnum} \
singularity exec --bind /hpcnfs/ \
${mergepeaks} \
${norm_compress_entropy_script}

#####idr
SINGULARITYENV_outdir=${outdir} SINGULARITYENV_rep1bam=${rep1bam} SINGULARITYENV_rep1bed=${rep1bed} \
SINGULARITYENV_rep2bam=${rep2bam} SINGULARITYENV_rep2bed=${rep2bed} SINGULARITYENV_inp1bam=${inp1bam} \
SINGULARITYENV_inp2bam=${inp2bam} SINGULARITYENV_inp1bed=${inp1bed} SINGULARITYENV_inp2bed=${inp2bed} \
SINGULARITYENV_inp1_readnum=${inp1_readnum} SINGULARITYENV_inp2_readnum=${inp2_readnum} \
SINGULARITYENV_rep1_readnum=${rep1_readnum} SINGULARITYENV_rep2_readnum=${rep2_readnum} \
singularity exec --bind /hpcnfs/ \
${idr} \
${idr_script}

###second norm
SINGULARITYENV_outdir=${outdir} SINGULARITYENV_rep1bam=${rep1bam} SINGULARITYENV_rep1bed=${rep1bed} \
SINGULARITYENV_rep2bam=${rep2bam} SINGULARITYENV_rep2bed=${rep2bed} SINGULARITYENV_inp1bam=${inp1bam} \
SINGULARITYENV_inp2bam=${inp2bam} SINGULARITYENV_inp1bed=${inp1bed} SINGULARITYENV_inp2bed=${inp2bed} \
SINGULARITYENV_inp1_readnum=${inp1_readnum} SINGULARITYENV_inp2_readnum=${inp2_readnum} \
SINGULARITYENV_rep1_readnum=${rep1_readnum} SINGULARITYENV_rep2_readnum=${rep2_readnum} \
singularity exec --bind /hpcnfs/ \
${mergepeaks} \
${norm_using_idr_peaks_script}















