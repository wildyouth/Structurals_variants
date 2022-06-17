#!/bin/bash

path_ref_gen='/users/gev/kot/stage/mauve/genome/Fasta_genome/genome_assemblies_genome_gb_PRJNA13758_N2_2013_fasta/GCA_000002985.3_WBcel235_genomic.fna'

if [[ ! -f ${path_ref_gen}.fai ]]
then
	samtools faidx $path_ref_gen
	cut -f1,2 ${path_ref_gen}.fai > ${path_ref_gen}_chrom_size.txt
fi

python3 /users/gev/kot/stage/mauve/csv_mauve_to_sv_compare.py $path_ref_gen 
