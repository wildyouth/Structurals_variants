#!/bin/bash

path='/users/gev/kot/stage/manta/manta_run_all/*'
path_end='/results/variants/diploidSV_inv.vcf.gz'
path_res='/users/gev/kot/stage/manta/all_bam_res/'

regions_file='/users/gev/kot/stage/manta/all_bam_res/regions.txt'

#for folder in $path;

#do echo ${folder}$path_end >> /users/gev/kot/stage/manta/run_all_res/diploidSV_vcf_list.txt;
#done;

cat $regions_file| while read line 
do
	#echo ${path_res}${line}_sv.txt $line

	python3 /users/gev/kot/stage/manta/svimmer-0.1/svimmer ${path_res}${line}_sv_inv.txt $line > ${path_res}${line}_inv.vcf
	
	bgzip ${path_res}${line}_inv.vcf
	tabix -p vcf ${path_res}${line}_inv.vcf.gz
done

#python3 /users/gev/kot/stage/manta/svimmer-0.1/svimmer ${path_res}diploidSV_vcf_list.txt CHROMOSOME_I CHROMOSOME_II CHROMOSOME_III CHROMOSOME_IV CHROMOSOME_V CHROMOSOME_X > ${path_res}all_bam_real.vcf



