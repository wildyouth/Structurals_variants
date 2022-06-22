#!/bin/bash

#script to perfom progressivemauve have the data store into csv file 

path_ref_gen=$1

path_target_gen=$2

path_output=$3

start_sw=$4

end_sw=$5

#Create a file with chromosomes size of reference genome in fasta format
if [[ ! -f ${path_ref_gen}.fai ]]
then
	samtools faidx $path_ref_gen
	cut -f1,2 ${path_ref_gen}.fai > ${path_ref_gen}_chrom_size.txt
fi

#run progressivemauve with sw_test.py
python3 /users/gev/kot/stage/mauve/sw_test.py $start_sw $end_sw $path_ref_gen $path_target_gen $path_output

#extract LCBs coordinate with their gaps to json file with mauve-parse.js
for xmfa in ${path_output}*.xmfa;
do node /users/gev/kot/stage/mauve/p3_mauve/scripts/mauve-parser.js --input $xmfa > $xmfa.json;
done;


#extract LCBs coordinate from json file to csv file and put a tag INV=inversion, INS=Insertion, INV_INS=Insertion and inversion for LCBs that are a SV
for json in  ${path_output}*.json;

do python3 /users/gev/kot/stage/mauve/json_to_csv.py $json;

done;

#extract only SV 
for csv in  ${path_output}/*.csv;

do python3 /users/gev/kot/stage/mauve/csv_mauve_to_sv_compare.py $csv ${path_ref_gen}_chrom_size.txt;

done;
