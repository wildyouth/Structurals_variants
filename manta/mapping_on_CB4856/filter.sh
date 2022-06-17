#!/bin/bash

bam_file = '/users/gev/kot/stage/manta/N2_and_CB4856_on_CB4856/mapping_on_CB4856/CB4856_on_CB4856_sorted.bam'

samtools view -h $bam_file | awk 'substr($0,1,1)=="@" || ($9>= 400)' |
samtools view -b > CB4856_on_CB4856_sorted_400_higher.bam
