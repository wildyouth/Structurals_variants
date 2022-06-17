#!/bin/bash
#bash script to run progressiveMauve
#Modify this file to put your files in the right folder

path='/users/gev/kot/stage/mauve/repeat/seed_weight_N2_WS220_vs_CB4856_19'
path_ref_gen='/users/gev/kot/stage/manta/c_elegans.WS220.dna.fa'

path_target_gen='/users/gev/kot/stage/mauve/genome/Fasta_genome/genome_assemblies_genome_gb_prjna523481_CB4856_2019_fasta/GCA_004526295.1_ASM452629v1_genomic.fna '

if [[ ! -f ${path_ref_gen}.fai ]]
then
	samtools faidx $path_ref_gen
	cut -f1,2 ${path_ref_gen}.fai > ${path_ref_gen}_chrom_size.txt
fi

for i in {18..24};
do
folder_name=${path}_repeat_$i;
mkdir $folder_name;

python3 /users/gev/kot/stage/mauve/sw_test.py 16 32 $path_ref_gen $path_target_gen ${folder_name}/


for xmfa in ${folder_name}/*.xmfa;
do node /users/gev/kot/stage/mauve/p3_mauve/scripts/mauve-parser.js --input $xmfa > $xmfa.json;
done;


for json in  ${folder_name}/*.json;
do python3 /users/gev/kot/stage/mauve/json_to_csv.py $json;
done;
mkdir ${folder_name}/csv;
mv ${folder_name}/*.csv;
for csv in  ${folder_name}/*.csv;
do python3 /users/gev/kot/stage/mauve/csv_mauve_to_sv_compare.py $path_ref_gen $csv;
done;
rm ${folder_name}/*.xmfa;
done;
