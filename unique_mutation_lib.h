#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//maximum sample limit, change it if you have more samples
#define MAXSAMPLE 1024
//maximum length of an indel
#define MAX_INDEL_LEN 128 
    
//bases are stored in arrays, element indices are defined here
#define ABASE 0
#define CBASE 1
#define GBASE 2
#define TBASE 3
#define REFBASE 4
#define DELBASE 5 

// 0 coverage samples also have some kind of "frequency"
// should always check for this when working with the frequencies
#define ZERO_COV_FREQ -9999

    
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
    int base_counts[MAXSAMPLE][6];
    int ins_counts[MAXSAMPLE];
    int del_counts[MAXSAMPLE];
   
    char** ins_bases[MAXSAMPLE];
    char** del_bases[MAXSAMPLE];
    
    double base_freqs[MAXSAMPLE][6];
    //these data are after base_Q filter, so the coverage can change
    int filtered_cov[MAXSAMPLE];
    
    //mutation info
    char mut_base;
    int mut_sample_idx;
    
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
    printf indels in mpileup line
*/
int print_mpileup_line_indels(struct Mpileup_line* my_pup_line);
    
/*
    printf base freqs in mpileup line
*/
int print_mpileup_line_freqs(struct Mpileup_line* my_pup_line);

/*
    print pos
*/
int print_mpileup_line_pos(struct Mpileup_line* my_pup_line);



////////////////////////////////////////////////////////////////////////////
// deep copy pileup struct
////////////////////////////////////////////////////////////////////////////

/*
    deep copy mpileup line
*/
int copy_mpileup_line(struct Mpileup_line* target ,struct Mpileup_line* from);


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
int count_bases_all_samples(struct Mpileup_line* my_pup_line,int baseq_lim);


/*
    counts bases in one sample
*/
int count_bases(char* bases, char* quals,
                int* base_counts,int* del_count, int* ins_count, int* filtered_cov,
                char ref_base,int baseq_lim);

/* 
    parse a base from the bases and quals
*/
int handle_base(char* bases,char* quals,int* base_counts, int* filtered_cov,
                   int* base_ptr,int* qual_ptr,int baseq_lim);


/* 
    parse a deletion from the bases and quals
*/
int handle_deletion(char* bases,int* del_count,
                   int* base_ptr);

/* 
    parse an insertion from the bases and quals
*/
int handle_insertion(char* bases,int* ins_count,
                   int* base_ptr);


////////////////////////////////////////////////////////////////////////////
// Collect indels
////////////////////////////////////////////////////////////////////////////

/*
    collect indels in all samples
*/
int collect_indels_all_samples(struct Mpileup_line* my_pup_line);


/*
    collect the inserted, and deleted bases
*/
int collect_indels(char* bases, char*** ins_bases, int ins_count, char*** del_bases, int del_count);

/*
    free memory of indel bases
*/
int free_indel_bases(char*** ins_bases,char*** del_bases);



////////////////////////////////////////////////////////////////////////////
// Proximal gap filtering related
////////////////////////////////////////////////////////////////////////////
  
/*
    updates last indel gap position if there is one at the position
*/
int update_last_gap(struct Mpileup_line* my_pup_line, char** last_gap_chrom, int* last_gap_pos, int* is_gap);

/*
    if position is gap delete all potential mutations too close
*/
int proximal_gap_hindsight_filter(struct Mpileup_line* potential_mut_lines,int* mut_ptr,
                                  char* last_gap_chrom,int last_gap_pos,
                                 int proximal_gap_min_distance);
    
/*
    prints and deletes mutations, which have survived the hindsight proximal gap filtering
*/
int flush_accepted_mutations(struct Mpileup_line* potential_mut_lines,
                             char* recent_chrom,int recent_pos,int* mut_ptr,
                             int proximal_gap_min_distance);

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
    print mutation
*/
int print_mutation(struct Mpileup_line* my_pup_line);

/*
    calls mutation from frequencies
*/
int call_mutation(struct Mpileup_line* my_pup_line,double sample_mut_freq_limit,
                  double min_other_ref_freq_limit,int cov_limit,
                  char* last_gap_chrom, int last_gap_pos, int proximal_gap_min_distance);

/*
    gets the highest not reference mut freq
*/
int get_max_non_ref_freq(struct Mpileup_line* my_pup_line,double* val, int* idx, char* mut_base);


/*
    gets the lowest reference freq, except for 1 sample
*/
int get_min_ref_freq(struct Mpileup_line* my_pup_line,double* val, int idx_2skip );

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