#!/bin/bash

file_list_read='/users/gev/kot/stage/manta/CB4856_manta/reads_inversions_chr1_43541_pos_read_name.txt'
bam_file='/users/gev/kot/stage/manta/N2_and_CB4856_on_CB4856/mapping_on_CB4856/CB4856_on_CB4856_sorted_400_higher.bam'

file=$(cat $file_list_read)

for line in $file
do
samtools view -h $bam_file | awk '$1=="$line" {print}'
echo $line
done

