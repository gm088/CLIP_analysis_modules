#!/bin/bash

WD=$PWD


#process for trimming and mapping of each sample
#replaces XXXX with the sample name/id thing in the template script
#replaces YYYY with a label of your choice

samples=()
labels=()
rawdatadir=""

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







