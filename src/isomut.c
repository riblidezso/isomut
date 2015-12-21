#include "isomut_lib.h"
    
#define MUT_BUFFER_SIZE 128  //must be longer than the proximal gap filtering distance! 
    
int main(int argc, char** argv)
{
    //cmdline args
    //parameters for mutation calling
    if(argc!=7){
        printf("ERROR please provide 5 args \n min_sample_freq (0.3) \n min_other_ref_freq (0.1)\n"); 
        printf(" cov_limit (10)\n base quality limit (30),\n prox_gap_dist_SNV  (10),\n  prox_gap_dist_indel  (100), \n"); 
        exit(1);
    }
    double min_sample_freq=strtod(argv[1],NULL);
    double max_other_freq=strtod(argv[2],NULL);
    int cov_limit=(int) strtol(argv[3],NULL,10);
    int baseq_limit=(int) strtol(argv[4],NULL,10);
    int prox_gap_min_dist_SNV= (int) strtol(argv[5],NULL,10);
    int prox_gap_min_dist_indel= (int) strtol(argv[6],NULL,10);
    
    //varaiables for reading a line
    char* line = NULL;
    size_t len = 0;
    ssize_t line_size;
    
    //potential mutation list
    struct mplp* potential_mutations = (struct mplp*) malloc(MUT_BUFFER_SIZE * sizeof(struct mplp)) ;
    //init all
    int i;
    for( i=0;i<MUT_BUFFER_SIZE;i++) init_mplp(&potential_mutations[i]);
    int mut_ptr = i = 0; //pointer for potential mutation buffer 
    
    //variables for proximal gap filtering
    int last_gap_pos_start,last_gap_pos_end;
    last_gap_pos_start = last_gap_pos_end = -42;
    char* last_gap_chrom = NULL;
    int is_gap = 1 ; // 0 yes, 1 no
    
    //print header
    printf("#sample_idx\tchr\tpos\ttype\tscore\tref\tmut\tcov\tmut_freq\tcleanliness\n");
    
    //loop over input lines
    //FILE* test_f = fopen("/Users/ribli/unique_mutation/data2/dbg.input","r");
    //while ((line_size = getline(&line, &len, test_f)) != -1) {
    while ((line_size = getline(&line, &len, stdin)) != -1) {
        //the pileup line structure for the line being read
        struct mplp my_mplp;
        init_mplp(&my_mplp);
        
        //build the struct from input line
        process_mplp_input_line(&my_mplp,line,line_size,baseq_limit);
        
        //call snvs with forward prox gap filtering
        call_snv(potential_mutations,&mut_ptr,&my_mplp,
                 min_sample_freq,max_other_freq,cov_limit,
                 last_gap_chrom,last_gap_pos_end,prox_gap_min_dist_SNV); 
        
        //call indels with forward prox gap filtering
        call_indel(potential_mutations,&mut_ptr,&my_mplp,
                   min_sample_freq,max_other_freq,cov_limit,
                   last_gap_chrom,last_gap_pos_start,last_gap_pos_end,
                   prox_gap_min_dist_SNV,prox_gap_min_dist_indel); 
                   
        //update last gap position seen
        update_last_gap(&my_mplp,&last_gap_chrom,&last_gap_pos_start,&last_gap_pos_end,&is_gap);
        
        //if a gap starts at this pos, delete in hindsight all potential mutations too close
        if( is_gap == 0 ) proximal_gap_hindsight_filter(potential_mutations,&mut_ptr,
                                                        last_gap_chrom,last_gap_pos_start,
                                                        prox_gap_min_dist_SNV,prox_gap_min_dist_indel);
        
        //if potential mut container is almost full: flush and print the accepted mutations
        if( mut_ptr > MUT_BUFFER_SIZE - 2 ){
            flush_accepted_mutations(potential_mutations,my_mplp.chrom,my_mplp.pos,
                                     &mut_ptr,prox_gap_min_dist_SNV,prox_gap_min_dist_indel);
        }
        
        //free the memory allocated by the struct
        free_mplp(&my_mplp);
    }
    
    //print mutations left at the very end
    for(i=0;i<mut_ptr;i++){ 
        print_mutation(&potential_mutations[i]);
        free_mplp(&potential_mutations[i]);
    }
    
    //free resources
    free(line);
    free(potential_mutations);
    return 0;
}