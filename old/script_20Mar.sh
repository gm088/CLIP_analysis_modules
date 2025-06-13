#!/bin/bash

pureclip="/Users/IEO5559/prog/clip/IP7_and_DCIP/pureclip/merged_reps/output_bed6.bed"
motif_bg="/Users/IEO5559/enhancedClip/promoter_bg_inverted_expr_14Mar.bed"
motif_bg2="/Users/IEO5559/enhancedClip/tx_5prime_nopeak_exprcontrol_randomsample_alt.bed"
motif_bg3="/Users/IEO5559/enhancedClip/MANE_promoters.bed"
motif_bg4="/Users/IEO5559/enhancedClip/MANE_promoters_inv.bed"
wd=`pwd`
pc_score=3


for file in `ls *_nohistone.bed`;
do
	#remove the ending things
	prefix=${file%_cpmfilt_gradfilt_nohistone.bed}
	#mkdir search_${prefix}
	cd search_${prefix}

	##remove empty peaks
	#intersectBed -v -s -a ../${file} -b ../empty_peaks_both.bed > ${prefix}_noempty.bed
#
#	###sort by score
#	#sort -k5,5nr ${prefix}_noempty.bed | head -n 5000 > top5000_byFDR.bed
#
#	###get rpkm scores
#	#Rscript ~/bin/get_rpkm_scores.R top5000_byFDR.bed yes 2> Rout	
#	##Rscript ~/bin/get_rpkm_scores.R ${prefix}_noempty.bed yes 2> Rout	
#
#	##sort -k5,5nr top5000_byFDR_rpkm.bed | head -n 1000 > top5000_byFDR_rpkm_top1000.bed
#	#sort -k5,5nr top5000_byFDR_rpkm.bed | head -n 500 > top5000_byFDR_rpkm_top500.bed
#
#	#rm ${prefix}_noempty.bed
#	##rm top5000_byFDR.bed
#
#	##Rscript ~/enhancedClip/peaksplitter_CLI.R top5000_byFDR_rpkm_top1000.bed min 2> Rout
#	##Rscript ~/enhancedClip/peaksplitter_CLI.R top5000_byFDR_rpkm_top1000.bed max 2> Rout
#
#	###overlap with PureCLIP XL sites 
	Rscript ~/bin/get_pureclip_xl.R ${pureclip} top5000_byFDR_rpkm_top500.bed ${pc_score} 2>> Rout

	##resulting file has XL score in column 4, and rpkm score in column 5

	##motif search using different backgrounds
	custom_bg_streme.sh 150 xl_bg1 \
	top5000_byFDR_rpkm_top500_xl_sites_collapsed5nt_${pc_score}.bed ${motif_bg} 4 6
	custom_bg_streme.sh 150 xl_bg2 \
	top5000_byFDR_rpkm_top500_xl_sites_collapsed5nt_${pc_score}.bed ${motif_bg2} 4 6
	custom_bg_streme.sh 150 xl_bg3 \
	top5000_byFDR_rpkm_top500_xl_sites_collapsed5nt_${pc_score}.bed ${motif_bg3} 4 6
	custom_bg_streme.sh 150 xl_bg4 \
	top5000_byFDR_rpkm_top500_xl_sites_collapsed5nt_${pc_score}.bed ${motif_bg4} 4 6

	#bottom xl bg
	sort -k5,5nr top5000_byFDR_rpkm.bed | tail -n 1000 > top5000_byFDR_rpkm_bottom1000.bed
	Rscript ~/bin/get_pureclip_xl.R ${pureclip} top5000_byFDR_rpkm_bottom1000.bed 0 2>> Rout
	bedtools slop -i top5000_byFDR_rpkm_bottom1000_xl_sites_collapsed5nt_0.bed \
	-g ~/human.hg38.genome -b 150 > bottomxl_bg.bed
	
	custom_bg_streme.sh 150 xl_xl_bg \
	top5000_byFDR_rpkm_top500_xl_sites_collapsed5nt_${pc_score}.bed \
	${wd}/search_${prefix}/bottomxl_bg.bed \
	4 6

	##distal bg
	Rscript ~/enhancedClip/bg_generator_distal.R \
	top5000_byFDR_rpkm_top500_xl_sites_collapsed5nt_${pc_score}.bed \
	1000 50 2>>Rout

	custom_bg_streme.sh 150 xl_distal_bg \
	top5000_byFDR_rpkm_top500_xl_sites_collapsed5nt_${pc_score}.bed \
	${wd}/search_${prefix}/top5000_byFDR_rpkm_top500_xl_sites_collapsed5nt_${pc_score}_distalbg.bed \
	4 6 


	##motif search in split peaks - top5000_byFDR_rpkm_top500_split_max.bed and min.bed
	###distal bg
	Rscript ~/enhancedClip/bg_generator_distal.R \
	top5000_byFDR_rpkm_top500_split_max.bed \
	1000 100 2>>Rout
	
	custom_bg_streme.sh 10 split_max_distal_bg \
	top5000_byFDR_rpkm_top500_split_max.bed \
	${wd}/search_${prefix}/top5000_byFDR_rpkm_top500_split_max_distalbg.bed \
	4 6 
	
	Rscript ~/enhancedClip/bg_generator_distal.R \
	top5000_byFDR_rpkm_top500_split_min.bed \
	1000 100 2>>Rout
	
	custom_bg_streme.sh 10 split_min_distal_bg \
	top5000_byFDR_rpkm_top500_split_min.bed \
	${wd}/search_${prefix}/top5000_byFDR_rpkm_top500_split_max_distalbg.bed \
	4 6 

	cd ..

done


















