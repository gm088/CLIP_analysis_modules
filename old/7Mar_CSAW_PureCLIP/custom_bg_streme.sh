#!/bin/bash

greed=$1
outdir=$2
bedfile=$3
bgbed=$4
min_width=$5
max_width=$6

#squash to TSS or TES
if [ -z $4 ]
then
  echo "usage: custom_bg_streme.sh greed outdir bedfile custom_bg_bed minw maxw"
  exit -1
fi

if [ -z $5 ]
then
  min_width=4
fi

if [ -z $6 ]
then
  max_width=6
fi

genome="/Users/IEO5559/Desktop/misc/splice/genome/genome.fa"
genome_size="/Users/IEO5559/human.hg38.genome"

mkdir $outdir 

cut -f 1-6 $bedfile > ${outdir}/peaks.bed

cd ${outdir}


##### CUSTOM BACKGROUND #######
bedtools getfasta -s -fi $genome -bed ${bgbed} -fo bg.fa

#eliminate sequences containing unmappable reds
awk '{a=$0; getline; if(!($1~ /N/)) {print a; print $0}}' bg.fa > bg_mappable.fa

#extract background for meme at different order
for m in 1 2 3 4 5; do
  if ! [ -a bg.order_${m}.param ]; then
    fasta-get-markov -m ${m} < bg_mappable.fa 1>bg.order_${m}.param 2>output.dat&
  fi
done

wait

#### slop
bedtools slop -i peaks.bed -g ~/human.hg38.genome -b ${greed} > tmppeaks.bed
##get sequences
bedtools getfasta -fi ${genome} -bed tmppeaks.bed -s > peakseqs.fa

#extract motifs with different background order
for o in 1 2 3 4; do
  
  streme --p peakseqs.fa --rna --minw $min_width --maxw $max_width --nmotifs 100 --order ${o} \
  -oc streme_${o} -bfile bg.order_${o}.param 2>nohup.out &
  
done

echo $1 $2 $3 $4 $min_width $max_width > params.txt





