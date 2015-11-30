#include "unique_mutation_lib.h"


////////////////////////////////////////////////////////////////////////////
// initialize
////////////////////////////////////////////////////////////////////////////
/*
    initializes pointers with NULL
*/
int init_mpileup_line(struct Mpileup_line* my_pup_line){
    int i;
    (*my_pup_line).raw_line=NULL;
    for(i=0;i<MAXSAMPLE;i++){
        (*my_pup_line).bases[i]=NULL;
        (*my_pup_line).quals[i]=NULL;
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////
// free resources
////////////////////////////////////////////////////////////////////////////
/*
    frees all malloced objects in Mpileup_line struct
*/
int free_mpileup_line(struct Mpileup_line* my_pup_line){
    free((*my_pup_line).raw_line);
    free((*my_pup_line).chrom);
    int i;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        free((*my_pup_line).bases[i]);
        free((*my_pup_line).quals[i]);
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////
// formatted printing functions
////////////////////////////////////////////////////////////////////////////


/*
    print the raw line
*/
int print_mpileup_raw_line(struct Mpileup_line* my_pup_line){
    //print the raw line
    printf("%s \n",(*my_pup_line).raw_line);
    return 0;
}



/*
    printf formatted representation of mpileup line
*/
int print_mpileup_line(struct Mpileup_line* my_pup_line){
    //print the raw line
    //printf("%s \n",(*my_pup_line).raw_line);

    //print position level info
    printf("%s %d %c\n",(*my_pup_line).chrom,(*my_pup_line).pos,(*my_pup_line).ref_nuq);
    
    //print bases and covs
    int i;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        printf("A %d,C %d,G %d,T %d,del %d,ins_start %d, del_start %d, cov %d\n",
               (*my_pup_line).base_counts[i][ABASE],
               (*my_pup_line).base_counts[i][CBASE],
               (*my_pup_line).base_counts[i][GBASE],
               (*my_pup_line).base_counts[i][TBASE],
               (*my_pup_line).base_counts[i][DELBASE],
               (*my_pup_line).ins_counts[i],
               (*my_pup_line).del_counts[i],
               (*my_pup_line).filtered_cov[i]);
                                        
        printf("A %.2f,C %.2f,G %.2f,T %.2f\n",(*my_pup_line).base_freqs[i][ABASE],
                                                (*my_pup_line).base_freqs[i][CBASE],
                                                (*my_pup_line).base_freqs[i][GBASE],
                                                (*my_pup_line).base_freqs[i][TBASE]);
        printf("%d %s %s\n",(*my_pup_line).cov[i],(*my_pup_line).bases[i],(*my_pup_line).quals[i]);
    }
    return 0;
}

/*
    printf base counts in mpileup line
*/
int print_mpileup_line_counts(struct Mpileup_line* my_pup_line){
    //print position level info
    printf("%s %d %c\n",(*my_pup_line).chrom,(*my_pup_line).pos,(*my_pup_line).ref_nuq);
    
    //print bases and covs
    int i;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        printf("A %d,C %d,G %d,T %d\n",(*my_pup_line).base_counts[i][ABASE],
                                        (*my_pup_line).base_counts[i][CBASE],
                                        (*my_pup_line).base_counts[i][GBASE],
                                        (*my_pup_line).base_counts[i][TBASE]);
    }
    return 0;
}


/*
    printf base freqs in mpileup line
*/
int print_mpileup_line_freqs(struct Mpileup_line* my_pup_line){
    //print position level info
    printf("%s %d %c\n",(*my_pup_line).chrom,(*my_pup_line).pos,(*my_pup_line).ref_nuq);
    
    //print bases and covs
    int i;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        printf("A %.2f,C %.2f,G %.2f,T %.2f\n",(*my_pup_line).base_freqs[i][ABASE],
                                                (*my_pup_line).base_freqs[i][CBASE],
                                                (*my_pup_line).base_freqs[i][GBASE],
                                                (*my_pup_line).base_freqs[i][TBASE]);
    }
    return 0;
}

/*
    print pos
*/
int print_mpileup_line_pos(struct Mpileup_line* my_pup_line){
    //print position level info
    printf("%s %d %c\n",(*my_pup_line).chrom,(*my_pup_line).pos,(*my_pup_line).ref_nuq);
    return 0;
}



////////////////////////////////////////////////////////////////////////////
// read pileup struct from mpileup line
////////////////////////////////////////////////////////////////////////////

/*
    gets mpileup line from the char* line
*/
int get_mpileup_line(struct Mpileup_line* my_pup_line,char* line, ssize_t line_size){   
    //store the raw line too
    free((*my_pup_line).raw_line);
    (*my_pup_line).raw_line = (char*)malloc( (line_size) * sizeof(char));
    memcpy((*my_pup_line).raw_line,line,(line_size) * sizeof(char));
    
    //temp buffer for reading those entries which will be formatted as not strings
    char* tmp_str=NULL;
    
    ssize_t i=0;
    while(i<line_size){
        //chrom
        get_next_entry(line,line_size,&i,&((*my_pup_line).chrom));
        //position
        get_next_entry(line,line_size,&i,&tmp_str);
        (*my_pup_line).pos=atoi(tmp_str);
        //ref nuq
        get_next_entry(line,line_size,&i,&tmp_str);
        (*my_pup_line).ref_nuq=tmp_str[0];
        
        //read samples
        int temp_sample=0;
        while(i<line_size){
            //coverage
            get_next_entry(line,line_size,&i,&tmp_str);
            (*my_pup_line).cov[temp_sample]=atoi(tmp_str);
            //bases
            get_next_entry(line,line_size,&i,&((*my_pup_line).bases[temp_sample]));
            //quals
            get_next_entry(line,line_size,&i,&((*my_pup_line).quals[temp_sample]));
            temp_sample++;
        }
        //store the number of samples
        (*my_pup_line).n_samples=temp_sample;
    }
    //free resources
    free(tmp_str);
    return 0;
}

/*
    gets next entry from tab separated input string as null terminated char*
*/
int get_next_entry(char* line, ssize_t line_size, ssize_t* pointer, char** result){
    int c,c0;
    c0 = *pointer;
    while(line[*pointer]!='\t' && *pointer<line_size)(*pointer)++;
    c = *pointer-c0;
    (*pointer)++;
    
    if(*result != NULL) free(*result);
    *result = (char*)malloc( (c+1) * sizeof(char));
    memcpy(*result,line+c0,c * sizeof(char));
    (*result)[c] = 0;
    
    //printf("%s \n",*result);
    return c;
    }



////////////////////////////////////////////////////////////////////////////
// Count bases in mplieup struct
////////////////////////////////////////////////////////////////////////////
/*
    counts bases in all samples
*/
int count_bases_all_samples(struct Mpileup_line* my_pup_line,int baseq_lim){
    int i=0;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        count_bases((*my_pup_line).bases[i],(*my_pup_line).quals[i],
                    (*my_pup_line).base_counts[i],
                    &((*my_pup_line).del_counts[i]),
                    &((*my_pup_line).ins_counts[i]),
                    &((*my_pup_line).filtered_cov[i]),
                    (*my_pup_line).ref_nuq,baseq_lim);
    }
    return 0;
}

/*
    counts bases in one sample
*/
int count_bases(char* bases, char* quals,int* base_counts, int* del_count,
                int* ins_count,int* filtered_cov,char ref_base,int baseq_lim){
    
    //initialize counts to zero
    base_counts[0] = base_counts[1] = base_counts[2] = 0;
    base_counts[3] = base_counts[4] = base_counts[5] = 0 ;
    *ins_count = *del_count = 0;
    
    int i=0; //pointer in str for bases
    int j=0; //pointer in str for qualities
    (*filtered_cov)=0;
    while(bases[i]!=0){
        //beginning and end of the read signs
        if(bases[i] == '$' ) i++; //next
        else if(bases[i] == '^' ) {i+=2;} //jump next character (mapq too)
        //deletions
        else if(bases[i]=='-' ) handle_deletion(bases,del_count,&i);
        //insetions
        else if(bases[i]=='+' ) handle_insertion(bases,ins_count,&i); 
        //real base data
        else handle_base(bases,quals,base_counts,filtered_cov,&i,&j,baseq_lim);
    }
    //add refbase to corresponding base
    if(ref_base=='A')base_counts[ABASE]+=base_counts[REFBASE];
    else if(ref_base=='C')base_counts[CBASE]+=base_counts[REFBASE];
    else if(ref_base=='G')base_counts[GBASE]+=base_counts[REFBASE];
    else if(ref_base=='T')base_counts[TBASE]+=base_counts[REFBASE];
    
    return 0;
}


/* 
    parse a base from the bases and quals
*/
int handle_base(char* bases,char* quals,int* base_counts, int* filtered_cov,
                   int* base_ptr,int* qual_ptr,int baseq_lim){
    if(quals[*qual_ptr] >= baseq_lim + 33 ){
        char c = bases[*base_ptr]; 
        if(c=='.' || c==',' )      base_counts[REFBASE]++;
        else if(c=='A' || c=='a' ) base_counts[ABASE]++;
        else if(c=='C' || c=='c' ) base_counts[CBASE]++;
        else if(c=='G' || c=='g' ) base_counts[GBASE]++;
        else if(c=='T' || c=='t' ) base_counts[TBASE]++;
        else if(c=='*' ) base_counts[DELBASE]++;
        (*filtered_cov)++;
    }
    (*qual_ptr)++;
    (*base_ptr)++;
    
    return 0;
}

/* 
    parse a deletion from the bases and quals
*/
int handle_deletion(char* bases,int* del_count,int* base_ptr){
    char* offset;
    (*del_count)++;
    int indel_len=strtol(&bases[*base_ptr+1],&offset,10);
    (*base_ptr)+= offset-&bases[*base_ptr+1] + indel_len;
 
    /*
    //save if first
    if( *del_count ==0 ) memcpy(del_bases,offset,indel_len * sizeof(char));
    //compare if second
    else{
        int i;
        for (i=0;i<indel_len;i++){
            if( offset[i] != del_bases[i] )
                
        }
     */   
    return 0;
}

/* 
    parse an insertion from the bases and quals
*/
int handle_insertion(char* bases,int* ins_count,int* base_ptr){
    char* offset;
    (*ins_count)++;
    int indel_len=strtol(&bases[*base_ptr+1],&offset,10);
    (*base_ptr)+= offset-&bases[*base_ptr+1] + indel_len;
    
    return 0;
}


////////////////////////////////////////////////////////////////////////////
// Calculate base frequencies
////////////////////////////////////////////////////////////////////////////

/*
    calculate base freqs in all samples
*/
int calculate_base_freqs_all_sample(struct Mpileup_line* my_pup_line){
    int i=0;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        calculate_base_freqs((*my_pup_line).base_freqs[i],
                            (*my_pup_line).base_counts[i],
                            (*my_pup_line).filtered_cov[i]);
    }
    return 0;
}


/*
   calculate base freqs in a sample
*/
int calculate_base_freqs(double* base_freqs,int* base_counts, int coverage){
    base_freqs[0] = base_freqs[1] = base_freqs[2] = base_freqs[3] = base_freqs[4] = 0 ;
    int i;
    if (coverage!=0){
        for(i=0;i<5;i++){
            base_freqs[i]=( (double) base_counts[i]) / coverage;
        }
    }
    else{
        for(i=0;i<5;i++){
            base_freqs[i]= ZERO_COV_FREQ;
        }
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////
// Call mutatations
////////////////////////////////////////////////////////////////////////////


/*
    calls mutation from frequencies
*/
int call_mutation(struct Mpileup_line* my_pup_line,double sample_mut_freq_limit,
                    double min_other_ref_freq_limit,int cov_limit){
    
    double sample_mut_freq,min_other_ref_freq;
    int sample_idx;
    char mut_base;
    
    //filter position if it has indel
    //todo
    
    get_max_non_ref_freq(my_pup_line,&sample_mut_freq,&sample_idx,&mut_base);
    get_min_ref_freq(my_pup_line,&min_other_ref_freq,sample_idx);

    if (sample_mut_freq >= sample_mut_freq_limit && // mut freq larger than limit
        min_other_ref_freq > min_other_ref_freq_limit && //other sample ref_freq higher than limit 
        (*my_pup_line).filtered_cov[sample_idx] >= cov_limit ){ //coverage higher than limit
        printf("%d\t%s\t%d\t%c\t%c\n",sample_idx,(*my_pup_line).chrom,(*my_pup_line).pos,
                (*my_pup_line).ref_nuq,mut_base);
    }
    return 0;        
}


/*
    gets the highest not reference mut freq
*/
int get_max_non_ref_freq(struct Mpileup_line* my_pup_line,double* val, int* idx, char* mut_base){

    //arrays only for encoding/decoding base indices to bases
    int base_2_idx[256];
    base_2_idx[ 'A' ] = 0;base_2_idx[ 'C' ] = 1;
    base_2_idx[ 'G' ] = 2;base_2_idx[ 'T' ] = 3;
    char idx_2_base[4];
    idx_2_base[ 0 ] = 'A';idx_2_base[ 1 ] = 'C';
    idx_2_base[ 2 ] = 'G';idx_2_base[ 3 ] = 'T';
  
    //get index for ref base
    int ref_idx=base_2_idx[(int)(*my_pup_line).ref_nuq];
   
    //initialize values to negative numbers
    *val = -42;
    *idx  = -42;
    
    int i,j;
    //loop over samples
    for(i=0;i<(*my_pup_line).n_samples;i++){
        //loop over bases
        for(j=0;j<4;j++){
            if ( j!=ref_idx && //base is not the reference base
                (*my_pup_line).base_freqs[i][j] > *val && //larger than largest yet
                (*my_pup_line).base_freqs[i][j] != ZERO_COV_FREQ ){ //not a 0 cov sample
                //save value of max, sample idx, and the mut base as chr
                *val=(*my_pup_line).base_freqs[i][j];
                *idx=i;
                *mut_base=idx_2_base[j];
            }
        }
    }
    return 0;
}


/*
    gets the lowest reference freq, except for 1 sample
*/
int get_min_ref_freq(struct Mpileup_line* my_pup_line,double* val, int idx_2skip ){
    //initialize value to large number
    *val = 42;
    
    int i;
    //loop over samples
    for(i=0;i<(*my_pup_line).n_samples;i++){
        if ( i != idx_2skip && // skip the mutated sample
           (*my_pup_line).base_freqs[i][REFBASE] < *val && //smaller than smallest yet
           (*my_pup_line).base_freqs[i][REFBASE] != ZERO_COV_FREQ ){ //not a 0 cov sample
                //save value of min
            *val=(*my_pup_line).base_freqs[i][REFBASE];
        }
    }
    return 0;
}

