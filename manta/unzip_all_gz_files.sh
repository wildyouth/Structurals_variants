#!/bin/bash

#unzip vcf.gz files

path='/users/gev/kot/stage/manta/test/manta_run_all/*'
path_end='/results/variants/diploidSV.vcf'


for folder in $path;

do gunzip -c ${folder}${path_end}.gz > ${folder}$path_end;


#bgzip ${folder};

#tabix -p vcf ${folder}.gz
done;

