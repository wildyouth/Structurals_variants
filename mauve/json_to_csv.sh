#!/bin/bash
#bash script to run progressiveMauve
#Modify this file to put your files in the right folder

path='/users/gev/kot/stage/mauve/seed_weed_ws220_cb4856/json'
path_ref_gen='/users/gev/kot/stage/manta/c_elegans.WS220.dna.fa'

path_target_gen='/users/gev/kot/stage/mauve/genome/Fasta_genome/genome_assemblies_genome_gb_prjna523481_CB4856_2019_fasta/GCA_004526295.1_ASM452629v1_genomic.fna '

folder_name='seed_weight_N2_WS220_vs_CB4856_19_repeat_'
#nb='1'
#for nb in {1..100};
#do
#for json in ${folder_name}${nb}/*.json;
for json in ${path}/*.json;
do python3 /users/gev/kot/stage/mauve/json_to_csv.py $json;
done;
#mv  ${folder_name}${nb}/*.csv;

for csv in  ${path}/*.csv;
do python3 /users/gev/kot/stage/mauve/csv_mauve_to_sv_compare.py $path_ref_gen $csv;
done;
#done;
