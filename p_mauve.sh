#!/bin/bash
#bash script to run progressiveMauve
#Modify this file to put your files in the right folder

path="/users/gev/kot/stage/mauve/Mauve_result/"
echo $path

name="N2_2013_vs_CB4856_2019" 

path_g1='/users/gev/kot/stage/mauve/genome/Fasta_genome/genome_assemblies_genome_gb_PRJNA13758_N2_2013_fasta/GCA_000002985.3_WBcel235_genomic.fna'
path_g2='/users/gev/kot/stage/mauve/genome/Fasta_genome/genome_assemblies_genome_gb_prjna523481_CB4856_2019_fasta/GCA_004526295.1_ASM452629v1_genomic.fna'
#path_g3='/users/gev/kot/Stage/Mauve/genome/GB_genome/genome_assemblies_genome_gb_prjna523481_CB4856_2019/GCA_004526295.1_ASM452629v1_genomic.gbk'
#path_g4='/users/gev/kot/Stage/Mauve/genome/GB_genome/genome_assemblies_genome_gb_CB4856_2015_gb/GCA_000975215.1_Cael_CB4856_1.0_genomic.gbk'
#path_g5='/users/gev/kot/Stage/Mauve/genome/GB_genome/genome_assemblies_genome_gb_PRJNA13758_N2_2014_gb/GCA_000939815.1_C_elegans_Bristol_N2_v1_5_4_genomic.gbk'
#path_g6 = 
#path_g7 = 
#path_g8 = 
#path_g9 = 

SECONDS=0

/export/home1/users/gev/kot/stage/mauve/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve --output=$path$name.xmfa --output-guide-tree=$path$name.tree --backbone-output=$path$name.backbone $path_g1 $path_g2

if (( $SECONDS > 3600 )) ; then
    let "hours=SECONDS/3600"
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo "Completed in $hours hour(s), $minutes minute(s) and $seconds second(s)" 
elif (( $SECONDS > 60 )) ; then
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo "Completed in $minutes minute(s) and $seconds second(s)"
else
    echo "Completed in $SECONDS seconds"
fi
'''
SECONDS=0

stripSubsetLCBs $path$name.xmfa $path$name.xmfa.bbcols $path$name\_stripSubsetLCBs.xmfa

if (( $SECONDS > 3600 )) ; then
    let "hours=SECONDS/3600"
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo "Completed in $hours hour(s), $minutes minute(s) and $seconds second(s)" 
elif (( $SECONDS > 60 )) ; then
    let "minutes=(SECONDS%3600)/60"
    let "seconds=(SECONDS%3600)%60"
    echo "Completed in $minutes minute(s) and $seconds second(s)"
else
    echo "Completed in $SECONDS seconds"
fi'''
