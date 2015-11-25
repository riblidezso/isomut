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
        printf("A %d,C %d,G %d,T %d, cov %d\n",(*my_pup_line).base_counts[i][ABASE],
                                        (*my_pup_line).base_counts[i][CBASE],
                                        (*my_pup_line).base_counts[i][GBASE],
                                        (*my_pup_line).base_counts[i][TBASE],
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
int count_bases_all_sample(struct Mpileup_line* my_pup_line,int baseq_lim){
    int i=0;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        count_bases((*my_pup_line).bases[i],(*my_pup_line).quals[i],
                    (*my_pup_line).base_counts[i],
                    &((*my_pup_line).filtered_cov[i]),
                    (*my_pup_line).ref_nuq,baseq_lim);
    }
    return 0;
}

/*
    counts bases in one sample
*/
int count_bases(char* bases,char* quals,int* base_counts,int* filtered_cov, char ref_base,int baseq_lim){
    base_counts[0] = base_counts[1] = base_counts[2] = base_counts[3] = base_counts[4] = 0 ;
    int i=0; //for bases
    int j=0; //for qualities
    char* offset;
    (*filtered_cov)=0;
    while(bases[i]!=0){
        //jump deletions insertions
        if(bases[i]=='-' || bases[i]=='+') {
            int indel_len=strtol(&bases[i+1],&offset,10);
            i+= offset-&bases[i+1] + indel_len;
            j--;
        }
        //jump the end of the read, beggining of the read signs
        else if(bases[i] == '^' ) {i++;j--;}
        else if(bases[i] == '$' ) j--;
       
        //real base data
        else if(quals[j] >= baseq_lim + 33 ){
            if(bases[i]=='.' || bases[i]==',' ) base_counts[REFBASE]++;
            else if(bases[i]=='A' || bases[i]=='a' ) base_counts[ABASE]++;
            else if(bases[i]=='C' || bases[i]=='c' ) base_counts[CBASE]++;
            else if(bases[i]=='G' || bases[i]=='g' ) base_counts[GBASE]++;
            else if(bases[i]=='T' || bases[i]=='t' ) base_counts[TBASE]++;
            (*filtered_cov)++;
            }
        i++;j++;
        }
        
    //add refbase to corresponding base
    if(ref_base=='A')base_counts[ABASE]+=base_counts[REFBASE];
    else if(ref_base=='C')base_counts[CBASE]+=base_counts[REFBASE];
    else if(ref_base=='G')base_counts[GBASE]+=base_counts[REFBASE];
    else if(ref_base=='T')base_counts[TBASE]+=base_counts[REFBASE];
    
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
            base_freqs[i]=-9999;
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

int call_mutation(struct Mpileup_line* my_pup_line,double min_sample_freq,
                    double max_other_freq,int cov_limit){
    double freq_1,freq_2;
    int idx_1,idx_2;
    char mut_base;
    get_2_highest_non_ref_freqs(my_pup_line,&freq_1,&freq_2,&idx_1,&idx_2,&mut_base);

    if (freq_1 >= min_sample_freq && freq_2 < max_other_freq && (*my_pup_line).filtered_cov[idx_1] >= cov_limit ){
        printf("%d\t%s\t%d\t%c\t%c\n",idx_1,(*my_pup_line).chrom,(*my_pup_line).pos,
                (*my_pup_line).ref_nuq,mut_base);
    }
    return 0;        
}


/*
    gets the 2 highest not reference base freqs
*/
int get_2_highest_non_ref_freqs(struct Mpileup_line* my_pup_line,double* val_1,
                                double* val_2, int* idx_1, int* idx_2, char* mut_base){
    if((*my_pup_line).n_samples<2){
        printf("ERROR get_2_highest_non_ref_freqs make no sense if there are less than 2 samples");
        exit(1);
    }
    int base_2_idx[256];
    base_2_idx[ 'A' ] = 0;base_2_idx[ 'C' ] = 1;
    base_2_idx[ 'G' ] = 2;base_2_idx[ 'T' ] = 3;
    char idx_2_base[4];
    idx_2_base[ 0 ] = 'A';idx_2_base[ 1 ] = 'C';
    idx_2_base[ 2 ] = 'G';idx_2_base[ 3 ] = 'T';
    
    int ref_idx=base_2_idx[(int)(*my_pup_line).ref_nuq];
    
    double tmp_val_1, tmp_val_2;
    tmp_val_1 = tmp_val_2 = -42;
    int tmp_idx_1, tmp_idx_2;
    tmp_idx_1 = tmp_idx_2 = -42;
    char tmp_mut_base=0;
    
    int i,j;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        for(j=0;j<4;j++){
            if ( j!=ref_idx && (*my_pup_line).base_freqs[i][j] > tmp_val_1 &&
                (*my_pup_line).base_freqs[i][j] !=-9999 ){
                tmp_val_1=(*my_pup_line).base_freqs[i][j];
                tmp_idx_1=i;
                tmp_mut_base=idx_2_base[j];
            }
        }
    }
    for(i=0;i<(*my_pup_line).n_samples;i++){
        for(j=0;j<4;j++){
            if ( j!=ref_idx && (*my_pup_line).base_freqs[i][j] > tmp_val_2 &&
                i!=tmp_idx_1 && (*my_pup_line).base_freqs[i][j] !=-9999 ){
                tmp_val_2=(*my_pup_line).base_freqs[i][j];
                tmp_idx_2=i;
            }
        }
    }
    
    *val_1=tmp_val_1;
    *val_2=tmp_val_2;
    *idx_1=tmp_idx_1;
    *idx_2=tmp_idx_2;
    *mut_base=tmp_mut_base;
    
    return 0;
}


/*
    gets the 2 smallest reference base freqs
    its not used now
*/
int get_2_smallest_ref_freqs(struct Mpileup_line* my_pup_line,double* val_1,
                                double* val_2, int* idx_1, int* idx_2){
    if((*my_pup_line).n_samples<2){
        printf("ERROR get_2_smallest_ref_freqs make no sense if there are less than 2 samples");
        exit(1);
    }
    
    double tmp_val_1, tmp_val_2;
    tmp_val_1 = tmp_val_2 = 42;
    int tmp_idx_1, tmp_idx_2;
    tmp_idx_1 = tmp_idx_2 = -42;
    
    int i;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        if ( (*my_pup_line).base_freqs[i][REFBASE] < tmp_val_1 &&
            (*my_pup_line).base_freqs[i][REFBASE] !=-9999 ){
            tmp_val_1=(*my_pup_line).base_freqs[i][REFBASE];
            tmp_idx_1=i;
        } 
    }
    for(i=0;i<(*my_pup_line).n_samples;i++){
        if ( (*my_pup_line).base_freqs[i][REFBASE] < tmp_val_2 &&
            (*my_pup_line).base_freqs[i][REFBASE] !=-9999 &&
            tmp_idx_1 !=i ){
            tmp_val_2=(*my_pup_line).base_freqs[i][REFBASE];
            tmp_idx_2=i;
        } 
    }
    
    *val_1=tmp_val_1;
    *val_2=tmp_val_2;
    *idx_1=tmp_idx_1;
    *idx_2=tmp_idx_2;
    return 0;
}



////////////////////////////////////////////////////////////////////////////
// Statistics
////////////////////////////////////////////////////////////////////////////

/*
    returns 1 if the position is clean in all samples , 0 if noisy
*/
int count_clean_pos(struct Mpileup_line* my_pup_line){
    int i;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        if((*my_pup_line).base_freqs[i][REFBASE] != 1.0 ) return 0;
    }
    return 1;
}
  

/*
    returns the number of sample covered with cov limit
*/
int count_covered_pos(struct Mpileup_line* my_pup_line, int cov_limit){
    int i,c;
    c=0;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        if((*my_pup_line).filtered_cov[i] >= cov_limit ) c++;
    }
    return c;
}

