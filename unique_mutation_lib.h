#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXSAMPLE 256

#define ABASE 0
#define CBASE 1
#define GBASE 2
#define TBASE 3
#define REFBASE 4

    
///////////////////////////////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////////////////////////////////    
    
/*
    the mpileup struct
*/
struct Mpileup_line
{
    //raw data
    char* raw_line;
    
    //position level data
    char* chrom;
    int pos;
    char ref_nuq;
    
    //sample data
    int n_samples;
    
    //sample level data
    int cov[MAXSAMPLE];
    char* bases[MAXSAMPLE];
    char* quals[MAXSAMPLE];
    
    //more structured data
    int base_counts[MAXSAMPLE][5];
    double base_freqs[MAXSAMPLE][5];
    //these data are after base_Q filter, so the coverage can change
    int filtered_cov[MAXSAMPLE];
    
};


////////////////////////////////////////////////////////////////////////////
// initialize
////////////////////////////////////////////////////////////////////////////
/*
    initializes pointers with NULL
*/
int init_mpileup_line(struct Mpileup_line* my_pup_line);



////////////////////////////////////////////////////////////////////////////
// free resources
////////////////////////////////////////////////////////////////////////////


/*
    frees all malloced objects in Mpileup_line struct
*/
int free_mpileup_line(struct Mpileup_line* my_pup_line);



////////////////////////////////////////////////////////////////////////////
// formatted printing functions
////////////////////////////////////////////////////////////////////////////


/*
    print the raw line
*/
int print_mpileup_raw_line(struct Mpileup_line* my_pup_line);


/*
    printf formatted representation of mpileup line
*/
int print_mpileup_line(struct Mpileup_line* my_pup_line);

/*
    printf base counts in mpileup line
*/
int print_mpileup_line_counts(struct Mpileup_line* my_pup_line);


/*
    printf base freqs in mpileup line
*/
int print_mpileup_line_freqs(struct Mpileup_line* my_pup_line);

/*
    print pos
*/
int print_mpileup_line_pos(struct Mpileup_line* my_pup_line);



////////////////////////////////////////////////////////////////////////////
// read pileup struct from mpileup line
////////////////////////////////////////////////////////////////////////////


/*
    gets next entry from tab separated input string as null terminated char*
*/
int get_next_entry(char* line, ssize_t line_size, ssize_t* pointer, char** result);

/*
    gets mpileup line from the char* line
*/
int get_mpileup_line(struct Mpileup_line* my_pup_line,char* line, ssize_t read);


////////////////////////////////////////////////////////////////////////////
// Count bases in mplieup struct
////////////////////////////////////////////////////////////////////////////

/*
    counts bases in all samples
*/
int count_bases_all_sample(struct Mpileup_line* my_pup_line,int baseq_lim);


/*
    counts bases in one sample
*/
int count_bases(char* bases,char* quals,int* base_counts,int* filtered_cov, char ref_base,int baseq_lim);


////////////////////////////////////////////////////////////////////////////
// Calculate base frequencies
////////////////////////////////////////////////////////////////////////////


/*
    calculate base freqs in all samples
*/
int calculate_base_freqs_all_sample(struct Mpileup_line* my_pup_line);


/*
   calculate base freqs in a sample
*/
int calculate_base_freqs(double* base_freqs,int* base_counts, int coverage);






////////////////////////////////////////////////////////////////////////////
// Call mutatations
////////////////////////////////////////////////////////////////////////////




/*
    gets the 2 smallest reference base freqs
*/
int get_2_smallest_ref_freqs(struct Mpileup_line* my_pup_line,double* val_1,
                                double* val_2, int* idx_1, int* idx_2);
                                

/*
    gets the 2 highest not reference base freqs
*/
int get_2_highest_non_ref_freqs(struct Mpileup_line* my_pup_line,double* val_1,
                                double* val_2, int* idx_1, int* idx_2, char* mut_base);

/*
    calls mutation from frequencies
*/

int call_mutation(struct Mpileup_line* my_pup_line,double min_sample_freq,
                    double max_other_freq,int cov_limit);




////////////////////////////////////////////////////////////////////////////
// Statistics
////////////////////////////////////////////////////////////////////////////


/*
    returns 1 if the position is clean in all samples , 0 if noisy
*/
int count_clean_pos(struct Mpileup_line* my_pup_line);

/*
    returns the number of sample covered with cov limit
*/
int count_covered_pos(struct Mpileup_line* my_pup_line, int cov_limit);