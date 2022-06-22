#!/bin/bash
#script to perfom progressivemauve multiple time and have the data store into csv file 

#path of reference genome
path_ref_gen=$1

#path of target genome
path_target_gen=$2

#path of repeat location to store data
path_repeat=$3

#staring seed_weight to test
start_sw=$4

#ending seed_weight to test
end_sw=$5

#number of repetition
nb_repeat=$6

#loop nb_repeat time
for ((i = 1; i <= $nb_repeat; i++));

do

#create folder to store data
folder_name=${path_repeat}repeat_${i}/;

mkdir $folder_name;

bash stage/mauve/json.sh $path_ref_gen $path_target_gen $folder_name $start_sw $end_sw;

#remove xmfa file because they are heavy and not needed after
rm ${folder_name}*.xmfa;

done;

