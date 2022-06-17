#!/bin/bash

ref_gen_fa='/users/gev/kot/stage/manta/c_elegans.WS220.dna.fa'

input_vcf='/users/gev/kot/stage/manta/all_bam_res/'

bam_list='/users/gev/kot/stage/manta/all_bam_res/'

regions_file='/users/gev/kot/stage/manta/all_bam_res/regions.txt'

output_path='/users/gev/kot/stage/manta/all_bam_res/results_gt/'

#do bgzip myvcf.vcf
#then tabix -p vcf myvcf.vcf.gz

cat $regions_file| while read line  
do
	/users/gev/kot/stage/manta/graphtyper genotype_sv $ref_gen_fa ${input_vcf}${line}_inv.vcf.gz --sams=${bam_list}${line}_bam_path.txt --region=$line --output=$output_path --threads=16
done
#for region in $region;
#do
#echo $b_list;
#done

#/users/gev/kot/stage/manta/graphtyper genotype_sv $ref_gen_fa $input_vcf --sams=$bam_list --region_file=$region --output=$output_path
