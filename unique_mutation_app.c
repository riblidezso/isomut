#include "unique_mutation_lib.h"
    
#define MUT_BUFFER_SIZE 128  //must be longer than the proximal gap filtering distance! 
    
int main(int argc, char** argv)
{
    //cmdline args
    //parameters for mutation calling
    if(argc!=6){
        printf("ERROR please provide 5 args \n min_sample_freq (0.3) \n min_other_ref_freq (0.1)\n"); 
        printf(" cov_limit (10)\n base quality limit (30),\n prox_gap_min_dist  (10) \n"); 
        exit(1);
    }
    double min_sample_freq=strtod(argv[1],NULL);
    double max_other_freq=strtod(argv[2],NULL);
    int cov_limit=strtol(argv[3],NULL,10);
    int baseq_limit=strtol(argv[4],NULL,10);
    int prox_gap_min_dist=strtol(argv[5],NULL,10);
    
    //varaiables for reading a line
    char* line = NULL;
    size_t len = 0;
    ssize_t chars_read;
    
    //the pileup structure for the line being read
    struct Mpileup_line my_pup_line;
    init_mpileup_line(&my_pup_line);
    
    //potential mutation list
    struct Mpileup_line* potential_mut_lines;
    potential_mut_lines = (struct Mpileup_line*) malloc(MUT_BUFFER_SIZE * sizeof(struct Mpileup_line)) ;
    //init all
    int i;
    for( i=0;i<MUT_BUFFER_SIZE;i++) init_mpileup_line(&(potential_mut_lines[i]));
    int mut_ptr = i = 0; //pointer for potential mutation buffer 
    
    //variables for proximal gap filtering
    int last_gap_pos_start,last_gap_pos_end;
    last_gap_pos_start = last_gap_pos_end = -42;
    char* last_gap_chrom = NULL;
    int is_gap = 1 ; // 0 yes, 1 no
    
    //print header
    printf("#sample\tchr\tpos\ttype\tref_nuq\tmut\n");
    
    //loop over input lines
    while ((chars_read = getline(&line, &len, stdin)) != -1) {
        //read mpileup 
        get_mpileup_line(&my_pup_line,line,chars_read);
        //count bases
        count_bases_all_samples(&my_pup_line,baseq_limit);
        //calculate freqs
        calculate_freqs_all_samples(&my_pup_line);
        
        //collect indels
        collect_indels_all_samples(&my_pup_line);
        //update last gap position seen
        update_last_gap(&my_pup_line,&last_gap_chrom,&last_gap_pos_start,&last_gap_pos_end,&is_gap);
        
        //if a gap start at this pos, delete in hindsight all potential mutations too close
        if( is_gap == 0 ) proximal_gap_hindsight_filter(potential_mut_lines,&mut_ptr,
                                                        last_gap_chrom,last_gap_pos_start,
                                                       prox_gap_min_dist);
        
        //call snvs with forward prox gap filtering
        call_snv(&(potential_mut_lines[mut_ptr]),&mut_ptr,&my_pup_line,
                      min_sample_freq,max_other_freq,cov_limit,
                      last_gap_chrom,last_gap_pos_end,prox_gap_min_dist); 
        
        //call indels with forward prox gap filtering
        call_indel(&(potential_mut_lines[mut_ptr]),&mut_ptr,&my_pup_line,
                      min_sample_freq,max_other_freq,cov_limit,
                      last_gap_chrom,last_gap_pos_start,last_gap_pos_end,prox_gap_min_dist); 
        
        //if potential mut container is full: flush and print the accepted mutations
        // also flush occasionally 
        if( mut_ptr == MUT_BUFFER_SIZE ){// || my_pup_line.pos %1000 == 0 ){
            flush_accepted_mutations(potential_mut_lines,my_pup_line.chrom,my_pup_line.pos,
                                     &mut_ptr,prox_gap_min_dist);
        }
    }
   
    //print mutations left at the very end
    for(i=0;i<mut_ptr;i++){ 
        print_mutation(&(potential_mut_lines[i]));
        free_mpileup_line( &(potential_mut_lines[i]) );
    }

    //free resources
    free(line);
    return 0;
}