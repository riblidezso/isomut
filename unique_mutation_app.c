#include "unique_mutation_lib.h"
#include <math.h>
#include <math.h>
        
int main(int argc, char** argv)
{
    //cmdline args
    if(argc!=5){
        printf("ERROR please provide 4 args \n min_sample_freq (0.3) \n max_other_freq (0.1)\n"); 
        printf(" cov_limit (10)\n base quality limit (30) \n"); 
        exit(1);
    }
    //parameters for mutation calling
    double min_sample_freq=strtod(argv[1],NULL);
    double max_other_freq=strtod(argv[2],NULL);
    int cov_limit=strtol(argv[3],NULL,10);
    int baseq_limit=strtol(argv[4],NULL,10);
    
    //varaiables for reading a line
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    
    //the pileup structure
    struct Mpileup_line my_pup_line;
    init_mpileup_line(&my_pup_line);
    
    //print header
    printf("#sample\tchr\tpos\tref_nuq\tmut\n");
    
    //loop over input lines
    while ((read = getline(&line, &len, stdin)) != -1) {
        //read mpileup, count bases, calculate freq, and call mutation
        get_mpileup_line(&my_pup_line,line,read);
        count_bases_all_sample(&my_pup_line,baseq_limit);
        calculate_base_freqs_all_sample(&my_pup_line);
        call_mutation(&my_pup_line,min_sample_freq,max_other_freq,cov_limit);
        
    }

    //free resources
    free_mpileup_line(&my_pup_line);
    free(line);
    return 0;
}