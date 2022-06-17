#!/bin/bash

#!/bin/bash

#path_bam_files=$(wc -l < $1)
#n2=$(wc -l < $2)

path_bam_files="/mnt/data3/mallard/mapping_on_CB4856/*.bam"
path_ref_genome='/users/gev/kot/stage/mauve/genome/Fasta_genome/genome_assemblies_genome_gb_prjna523481_CB4856_2019_fasta/GCA_004526295.1_ASM452629v1_genomic.fna'
path_output='/users/gev/kot/stage/manta/N2_and_CB4856_on_CB4856/'
path_end='/results/variants/diploidSV'
#go through each bam file
for bam in $path_bam_files;

#setup folder name
do folder=$path_output$(basename "${bam%.*}");
#manta 
configManta.py --bam=$bam --referenceFasta=$path_ref_genome --runDir=$folder;


#analysis
$folder/runWorkflow.py;
#

gunzip -c ${folder}${path_end}.vcf.gz > ${folder}${path_end}.vcf


/users/gev/kot/miniconda3/envs/manta/share/manta-1.6.0-1/libexec/convertInversion.py /usr/bin/samtools $path_ref_genome ${folder}${path_end}.vcf > ${folder}${path_end}_inv.vcf

bgzip -c ${folder}${path_end}_inv.vcf > ${folder}${path_end}_inv.vcf.gz
tabix -p vcf ${folder}${path_end}_inv.vcf.gz
done;


#parallel -k "smartctl -a /dev/{}" ::: {/mnt/usb_ext/Noble_GDrive/bam_real/A*.bam} > path/to/output

