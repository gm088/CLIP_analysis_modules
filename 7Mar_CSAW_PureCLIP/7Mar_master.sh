#!/bin/bash

pureclip="/Users/IEO5559/prog/clip/IP7_and_DCIP/pureclip/merged_reps/output_bed6.bed"
motif_bg="/Users/IEO5559/enhancedClip/tx_5prime_nopeak_exprcontrol_randomsample.bed"
wd=`pwd`

for file in `ls *_nohistone.bed`;
do
	#remove the ending things
	prefix=${file%_cpmfilt_gradfilt_nohistone.bed}
	mkdir ${prefix} && cd ${prefix}

	##temporary fix
	cp ../${file} .

	##get rpkm scores
	Rscript ~/bin/get_rpkm_scores.R ${file} 2> Rout	

	##temporary fix part2
	rm ${file}

	##overlap with PureCLIP XL sites 
	Rscript ~/bin/get_pureclip_xl.R ${pureclip} ${file%.bed}_rpkm.bed 2>> Rout

	##resulting file has XL score in column 4, and rpkm score in column 5

	##for 3 different xl score thresh
	for i in {3..5};
	do
		#mkdir xl_${i} && cd xl_${i}
		#take xl score > 3,4,5
		awk -v var="$i" '$4>var' xl_sites.bed > xl_sites_gt${i}.bed	
		awk -v var="$i" '$4>var' xl_sites_collapsed1kb.bed > xl_sites_collapsed1kb_gt${i}.bed

		#sort by rpkm score
		sort -k5,5nr xl_sites_gt${i}.bed >  xl_sites_gt${i}_sorted.bed
		sort -k5,5nr xl_sites_collapsed1kb_gt${i}.bed >  xl_sites_collapsed1kb_gt${i}_sorted.bed
		rm  xl_sites_gt${i}.bed xl_sites_collapsed1kb_gt${i}.bed

		#### MESS BEGINS HERE

		head -n 500 xl_sites_gt${i}_sorted.bed > xl_sites_gt${i}_sorted_top500.bed
		head -n 500 xl_sites_collapsed1kb_gt${i}_sorted.bed > xl_sites_collapsed1kb_gt${i}_sorted_top500.bed
		##motif search
		custom_bg_streme.sh 300 xl_${i}_top500 xl_sites_gt${i}_sorted_top500.bed ${motif_bg}
		wait
		custom_bg_streme.sh 300 xl_${i}_top500_collapsed xl_sites_collapsed1kb_gt${i}_sorted_top500.bed ${motif_bg}
		wait

		head -n 1000 xl_sites_gt${i}_sorted.bed > xl_sites_gt${i}_sorted_top1000.bed
		head -n 1000 xl_sites_collapsed1kb_gt${i}_sorted.bed > xl_sites_collapsed1kb_gt${i}_sorted_top1000.bed 
		###motif search
		custom_bg_streme.sh 300 xl_${i}_top1000 xl_sites_gt${i}_sorted_top1000.bed ${motif_bg}
		wait
		custom_bg_streme.sh 300 xl_${i}_top1000_collapsed xl_sites_collapsed1kb_gt${i}_sorted_top1000.bed ${motif_bg}
		wait

	done

	cd ..

done


