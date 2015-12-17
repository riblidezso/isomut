from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing,sys,time
import subprocess
import glob

#
# Calculates the blocks of parallelization on the genome
#     - No blocks overlapping chromosomes
#          - gives back more than min_block_no blocks!
#
def define_parallel_blocks(ref_genome,min_block_no):
    print 'Defining parallell blocks ...'
    sys.stdout.flush()
    #genome lengths
    chromDict=dict()
    fullLeng=0

    #parsing fasta file
    for seqin in SeqIO.parse(ref_genome,"fasta"):
        # scaffolds and mitochondrium not used
        if (len(seqin.id) < 3 and seqin.id!='MT' ):
            #save lengths
            chromDict[seqin.id]=len(seqin.seq)
            fullLeng+=len(seqin.seq)

    #set maximum block size
    BLOCKSIZE=(fullLeng/min_block_no)

    #calculate blocks
    blocks=[]
    #loop pover chroms
    for chrom in sorted(chromDict.keys()):
        pointer=0
        #until chrom is chopped into pieces
        while (pointer < chromDict[chrom]):
            blockSize=min(BLOCKSIZE,chromDict[chrom]-pointer)
            blocks.append([chrom,pointer,pointer+blockSize])
            pointer += blockSize
    
    #reurn chr, posmin, posmax of blocks
    print 'Done'
    print
    sys.stdout.flush()
    return blocks



def run_isomut_on_block(chrom,from_pos,to_pos,
                         input_dir,bam_files,output_dir,
                         ref_genome,
                         min_sample_freq=0.2,
                         min_other_ref_freq=0.8,
                         cov_limit=10,
                         base_quality_limit=30,
                         min_gap_dist_snv=10,
                         min_gap_dist_indel=20,
                         bedfile=None,
                         samtools_flags=' -B '):
    #build the command
    cmd=' samtools  mpileup ' + samtools_flags
    cmd+=' -f ' +ref_genome 
    cmd+=' -r '+chrom+':'+str(from_pos)+'-'+str(to_pos)+' '
    if(bedfile!=None):
        cmd+=' -l '+bedfile+' '
    for bam_file in bam_files:
        cmd+=input_dir+bam_file +' '    
    cmd+=' 2> samtools.log | ./isomut '
    cmd+=' '.join(map(str,[min_sample_freq,min_other_ref_freq,cov_limit,
                           base_quality_limit,min_gap_dist_snv,min_gap_dist_indel])) +' '
    cmd+=' > ' +output_dir+'/'+ (chrom)+'_'+str(from_pos)+'_'+str(to_pos)+'_mut.csv  '
    
    return subprocess.check_call(cmd,shell=True)


def run_isomut_in_parallel(params):
    #check for bedfile argument:
    if (not params.has_key('bedfile')):
        params['bedfile']=None
    if (not params.has_key('samtools_flags')):
        params['samtools_flags']=' -B '
        
    #define blocks and create args
    blocks=define_parallel_blocks(params['ref_fasta'],params['n_min_block'])
    args=[]
    for block in blocks: 
        args.append([ block[0],block[1],block[2],
                     params['input_dir'],params['bam_filenames'],
                     params['output_dir'],params['ref_fasta'],
                     params['min_sample_freq'],params['min_other_ref_freq'],
                     params['cov_limit'], params['base_quality_limit'],
                     params['min_gap_dist_snv'],params['min_gap_dist_indel'],
                     params['bedfile'],params['samtools_flags']])
        
    #create dir
    if(glob.glob(params['output_dir'])==[]):
        subprocess.call(['mkdir',params['output_dir']])
    
    print 'blocks to run:',len(args)
    print 'running:',
    #start first n concurrent block
    procs=[]
    for i in xrange(params['n_conc_blocks']):
        procs.append(multiprocessing.Process(target=run_isomut_on_block, args=args[len(procs)]))
        procs[-1].start()
        print len(procs),
        sys.stdout.flush()
        
    # when one finished restart teh next
    while (len(procs) != len(args)):
        for i in xrange(1,params['n_conc_blocks']+1):
            #one finished start another one
            if(procs[-i].is_alive() == False and len(procs) != len(args)):
                procs.append(multiprocessing.Process(target=run_isomut_on_block,args=args[len(procs)] ))
                procs[-1].start()
                print len(procs),
        time.sleep(0.1)

    #wait now only the last ones running
    finished = False
    while( not finished):
        finished = True
        for proc in procs:
            if (proc.is_alive() == True):
                finished = False
        time.sleep(0.1)
        
    print '\nDone'
    
    
def run_isomut_with_pp(params):
    #run first
    run_isomut_in_parallel(params)

    # collect indels
    with open(params['output_dir']+'/all_indels.isomut','w') as indel_f  : 
        #write header
        indel_f.write("#sample\tchr\tpos\ttype\tscore\tref\tmut\tdepth\tmut_freq\n")
    #copy all indel lines except the header and sort them by chr, pos
    subprocess.check_call('cat ' +params['output_dir']+'/*.csv | grep -v "#"  | \
    grep -v SNV  | sort -n -k2,2 -k3,3 >> '+params['output_dir']+'/all_indels.isomut',shell=True)

    # create bedfile for SNVs for post processing
    subprocess.check_call('cat ' +params['output_dir']+'/*.csv | grep -v "#"  | \
    grep SNV | sort -n -k2,2 -k3,3 | cut -f 2,3 >'+params['output_dir']+'/temp.bed',shell=True)

    # clean everything else
    subprocess.check_call('rm '+params['output_dir']+'/*.csv',shell=True)

    # change params for postprocessing
    params['base_quality_limit']= 13
    params['min_other_ref_freq']= 0
    params['samtools_flags'] = ' ' 
    params['bedfile']=params['output_dir']+'/temp.bed'
    #run it
    run_isomut_in_parallel(params)


    # collect SNVs
    with open(params['output_dir']+'/all_SNVs.isomut','w') as snv_f:
        snv_f.write("#sample\tchr\tpos\ttype\tscore\tref\tmut\tdepth\tmut_freq\n")
    subprocess.check_call('cat ' +params['output_dir']+'/*.csv | grep -v "#"  | \
    grep SNV | sort -n -k2,2 -k3,3 >> '+params['output_dir']+'/all_SNVs.isomut',shell=True)

    #clean up
    subprocess.check_call(['rm',params['bedfile']])
    subprocess.check_call('rm '+params['output_dir']+'/*.csv',shell=True)
    
    print '\nDone'