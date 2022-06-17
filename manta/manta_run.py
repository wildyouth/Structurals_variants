# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import subprocess
import os
import glob
# If your shell script has shebang, 
# you can omit shell=True argument.

ref_fasta = '/users/gev/kot/stage/manta/c_elegans.WS220.dna.fa'

os.chdir(r'/users/gev/kot/stage/manta/bam_files/')
bam_files = glob.glob('*.bam')

for bam_f in bam_files:
    
    bam_n2 = '/users/gev/kot/stage/manta/bam_files/N2_ANC.bam'
    bam_cb4856 = '/users/gev/kot/stage/manta/bam_files/CB4856.A6330.bam'
    
    
    path_output_n2 = '/users/gev/kot/stage/manta/results/' + 'res_n2'

    argument = ' --bam ' + bam_n2 + ' --bam ' + bam_f + ' -- referenceFasta ' + ref_fasta + ' --runDir ' + path_output_n2

    subprocess.run(['python /users/gev/kot/stage/manta/manta-1.6.0.release_src/src/python/bin/configManta.py'])# + argument])
    
    ################################################"
    path_output = '/users/gev/kot/stage/manta/results/' + 'res_cb4856'

    argument = ' --bam ' + bam_n2 + ' --bam ' + bam_cb4856 + ' -- referenceFasta ' + ref_fasta + ' --runDir ' + path_output

    subprocess.run(['/bin/bash/python /users/gev/kot/stage/manta/manta-1.6.0.release_src/src/python/bin/configManta.py'])# + argument])