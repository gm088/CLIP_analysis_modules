#!/bin/bash

WD=$PWD

#make necessary directories

#mkdir FASTQ metrics trims mapped bigwigs subscripts bigwigs_norm repeat_mapped umi
#mkdir trims/qc

#process for trimming and mapping of each sample
#replaces XXXX with the sample name/id thing in the template script
#replaces YYYY with a label of your choice

samples=('Sample_S60783_PNUTS_con_WB_NT' 'Sample_S60784_PNUTS_con_WB_AUX' 'Sample_S61129_20250304_PNUTS_WT_R1' 'Sample_S61130_20250304_PNUTS_WT_R2' 'Sample_S61131_20250304_PNUTS_delRNA_R1' 'Sample_S61132_20250304_PNUTS_delRNA_R2' 'Sample_S61133_20250304_ZC3H4_WT_R4')
labels=('PNUTS_eCLIP_NT_6Feb25_2' 'PNUTS_eCLIP_AUX_6Feb25_2' 'PNUTS_delRNA_NT_R1' 'PNUTS_delRNA_NT_R2' 'PNUTS_delRNA_DOX_R1' 'PNUTS_delRNA_DOX_R2' 'WT_22May_R4_2')
#rawdatadir="/hpcnfs/techunits/genomics/PublicData/Natoli/dpolizzese/FASTQ/250204_A00302_0694_AH2T7JDMX2"
rawdatadir="/hpcnfs/techunits/genomics/PublicData/Natoli/dpolizzese/FASTQ/250317_A00302_0700_BHN5JGDRX5"

script='clip_libtype1.sh'

if [ ${#samples[@]} -ne ${#labels[@]} ]; then
  echo "samples and labels of unequal length"
  exit
fi

name=${script%.sh}

i=0

while [ ${i} -lt ${#samples[@]} ]
do
  
  numjobs=`qstat -u ieo5559 | grep "^[0-9]" | wc -l`
  if [ $numjobs -lt 6 ]
  then

    awk -v var1="${samples[i]}" -v var2="${labels[i]}" -v var3=${rawdatadir} '{gsub(/XXXX/, var1); gsub(/YYYY/, var2); gsub(/X_RAWDATADIR_X/, var3); print}' ${script} \
    > subscripts/${name}_${labels[i]}.sh
    qsub subscripts/${name}_${labels[i]}.sh
    i=$(( $i + 1 ))
    
  else
    
    sleep 1h
    
  fi
done







