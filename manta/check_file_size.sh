#!/bin/bash

path='/users/gev/kot/stage/manta/manta_run_all/*'

path_error='/users/gev/kot/stage/manta/error_bam/failed_bam_manta/'

for folder in $path;
do
if [ -s ${folder}/workflow.error.log.txt ];
then
    echo $folder
    
    #move failed bam
    mv $folder $path_error
    
#else
 #   echo "Size of sample.txt is zero"
fi
done

#do
#echo ${folder}/workflow.error.log.txt;
#s=$(wc -c ${folder}/workflow.error.log.txt | awk '{print $1}')
#echo $s
#if [-s $s];
#then
#echo 'youpi'
#fi
#done

