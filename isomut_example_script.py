#!/usr/bin/env python

#################################################
# importing the wrapper
#################################################
import sys,os
#add path for isomut_wrappers.py
#	if not running it from the isomut directory
#	change os.getcwd for the path to it
sys.path.append(os.getcwd()+'/src')
#load the parallel wrapper function
from isomut_wrappers import run_isomut

#add path for isomut, if its in the path comment/delete this line
#	if not running it from the isomut directory
#	change os.getcwd for the path to it
os.environ["PATH"] += os.pathsep + os.getcwd() +'/src'


#################################################
# defining administrative parameters
#################################################
#using parameter dictionary, beacause there are awful lot of parameters
params=dict()
#minimum number of blocks to run
# usually there will be 10-20 more blocks
params['n_min_block']=200
#number of concurrent processes to run
params['n_conc_blocks']=4
#genome
params['ref_fasta']="/home/ribli/input/index/gallus/Gallus_gallus.Galgal4.74.dna.toplevel.fa"
#input dir output dir
params['input_dir']='/nagyvinyok/adat86/sotejedlik/ribli/dt40/test_bams/'
params['output_dir']='isomut_test_output/'
#the bam files used
params['bam_filenames']=['DS014.bam', 'DS051.bam', 'DS052.bam', 'DS053.bam', 'DS054.bam', 'DS055.bam',
         'DS056.bam', 'DS057.bam', 'DS058.bam', 'DS101.bam', 'DS102.bam', 'DS103.bam']

#limit chromosomes (for references with many scaffolds)
# just comment/delete this line if you want to analyze all contigs in the ref genome
params['chromosomes']=map(str,range(1,29))+ ['32','W','Z','MT']

#################################################
# defining mutation calling parameters
#    default values here ...
#################################################
params['min_sample_freq']=0.21
params['min_other_ref_freq']=0.93
params['cov_limit']=5
params['base_quality_limit']=30
params['min_gap_dist_snv']=0
params['min_gap_dist_indel']=20

#################################################
# and finally run it
#################################################
run_isomut(params)
