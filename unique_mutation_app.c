#include "unique_mutation_lib.h"
#include <math.h>
#include <math.h>
        
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
    ssize_t read;
    
    //the pileup structure
    struct Mpileup_line my_pup_line;
    init_mpileup_line(&my_pup_line);
    
    //variables for proximal gap filtering
    int last_gap_pos = -42;
    char* last_gap_chrom = NULL;
    
    //print header
    printf("#sample\tchr\tpos\tref_nuq\tmut\n");
    
    //loop over input lines
    while ((read = getline(&line, &len, stdin)) != -1) {
        //read mpileup 
        get_mpileup_line(&my_pup_line,line,read);
        
        //count bases
        count_bases_all_samples(&my_pup_line,baseq_limit);
        
        //collect indels
        collect_indels_all_samples(&my_pup_line);
        
        //update last gap position seen
        update_last_gap(&my_pup_line,&last_gap_chrom,&last_gap_pos);
        
        //calculate freq, and call mutation
        calculate_base_freqs_all_sample(&my_pup_line);
        call_mutation(&my_pup_line,min_sample_freq,max_other_freq,cov_limit,
                     last_gap_chrom,last_gap_pos,prox_gap_min_dist);
       
        //printf("%s %d\n",last_gap_chrom,last_gap_pos);
        
        //print_mpileup_line(&my_pup_line);
        if( my_pup_line.del_counts[0]!= 0 ){
            //print_mpileup_line(&my_pup_line);
            //print_mpileup_line_indels(&my_pup_line);
        }
    }

    //free resources
    free_mpileup_line(&my_pup_line);
    free(line);
    return 0;
}