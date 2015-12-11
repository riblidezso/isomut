#!/usr/bin/env python

#add path for isomput_parallel.py if its not here
#import sys
#sys.path.append('path for isomut_parallel.py')

#load the parallel wraooer function
from isomut_parallel import run_isomut_in_parallel

#minimum number of blocks to run
# usually there will be 10-20 more blocks
min_block_no=200

#number of concurrent processes to run
no_conc_blocks=4

#genome
ref_fasta="/home/ribli/input/index/gallus/Gallus_gallus.Galgal4.74.dna.toplevel.fa"

#input dir output dir
input_dir='/nagyvinyok/adat86/sotejedlik/ribli/dt40/test_bams/'
output_dir='parall_test_output/'

#the bam files used
bam_filenames=['DS014.bam', 'DS051.bam', 'DS052.bam', 'DS053.bam', 'DS054.bam', 'DS055.bam',
         'DS056.bam', 'DS057.bam', 'DS058.bam', 'DS101.bam', 'DS102.bam', 'DS103.bam']

#and run it
run_isomut_in_parallel(input_dir,bam_filenames,output_dir,ref_fasta,
                       min_block_no=min_block_no,n_conc_blocks=no_conc_blocks)