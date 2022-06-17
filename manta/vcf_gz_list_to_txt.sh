#!/bin/bash

path='/users/gev/kot/stage/manta/manta_run_all/*'
path_end='/results/variants/diploidSV_inv.vcf'
path_res='/users/gev/kot/stage/manta/all_bam_res/'

for folder in $path;

do 

#echo ${folder}$path_end.gz >> ${path_res}diploidSV_vcf_inv_list.txt;

#gunzip -c ${folder}${path_end}.gz > ${folder}$path_end;
tabix -p vcf ${folder}${path_end}.gz
done;
