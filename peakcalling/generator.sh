#!/bin/bash

WD=$PWD

if [ -z $6 ]; then echo "usage: generator.sh outdir rep1bam rep1bed rep2bam rep2bed name_of_exp"; exit; fi

outdir=$1

rep1bam=$2
rep1bed=$3

rep2bam=$4
rep2bed=$5

name=$6

script='/hpcnfs/data/GN2/gmandana/bin/clipper_idr_scripts/post_clipper.sh'

awk -v var1="${outdir}" -v var2="${rep1bam}" -v var3="${rep1bed}" -v var4="${rep2bam}" -v var5="${rep2bed}" \
'{gsub(/OUTDIR/, var1); gsub(/REP1BAM/, var2); gsub(/REP1BED/, var3); gsub(/REP2BAM/, var4); gsub(/REP2BED/, var5); print}' ${script} \
    > norm_IDR_${name}.sh


echo ${outdir} ${rep1bam} ${rep1bed} ${rep2bam} ${rep2bed} ${name} > ${name}.params
