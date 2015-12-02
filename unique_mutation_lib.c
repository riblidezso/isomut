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
    (*my_pup_line).chrom=NULL;
    for(i=0;i<MAXSAMPLE;i++){
        (*my_pup_line).bases[i]=NULL;
        (*my_pup_line).quals[i]=NULL;
        (*my_pup_line).ins_bases[i]=NULL;
        (*my_pup_line).del_bases[i]=NULL;
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
    int i,j;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        free((*my_pup_line).bases[i]);
        free((*my_pup_line).quals[i]);
        
        for(j=0;j<(*my_pup_line).ins_counts[i];j++){
            free((*my_pup_line).ins_bases[i][j]);
        }
        free((*my_pup_line).ins_bases[i]);
        
        for(j=0;j<(*my_pup_line).del_counts[i];j++){
            free((*my_pup_line).del_bases[i][j]);
        }
        free((*my_pup_line).del_bases[i]);
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
    printf indels in mpileup line
*/
int print_mpileup_line_indels(struct Mpileup_line* my_pup_line){
    //print position level info
    printf("%s %d %c\n",(*my_pup_line).chrom,(*my_pup_line).pos,(*my_pup_line).ref_nuq);
    
    //print bases and covs
    int i,j;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        printf("\n\ninsertions: ");
        for(j=0;j<(*my_pup_line).ins_counts[i];j++)
            printf("%s ",(*my_pup_line).ins_bases[i][j]);
        printf("\ndeletions: ");
        for(j=0;j<(*my_pup_line).del_counts[i];j++)
            printf("%s ",(*my_pup_line).del_bases[i][j]);
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
// deep copy pileup struct
////////////////////////////////////////////////////////////////////////////

/*
    deep copy mpileup line
*/
int copy_mpileup_line(struct Mpileup_line* target ,struct Mpileup_line* from) {
    
    (*target).raw_line = (char*) malloc((strlen((*from).raw_line)+1) * sizeof(char));
    strcpy( (*target).raw_line, (*from).raw_line);
    
    (*target).chrom = (char*) malloc((strlen((*from).chrom)+1) * sizeof(char));
    strcpy( (*target).chrom, (*from).chrom);

    (*target).pos=(*from).pos;
    (*target).ref_nuq=(*from).ref_nuq;
    (*target).n_samples=(*from).n_samples;
    
    (*target).mut_base=(*from).mut_base;
    (*target).mut_sample_idx=(*from).mut_sample_idx;
   
    int i,j;
    for(i=0;i<(*target).n_samples;i++){
        (*target).cov[i]=(*from).cov[i];
        
        (*target).bases[i] = (char*) malloc((strlen((*from).bases[i])+1) * sizeof(char));
        strcpy( (*target).bases[i], (*from).bases[i]);
        
        (*target).quals[i] = (char*) malloc((strlen((*from).quals[i])+1) * sizeof(char));
        strcpy( (*target).quals[i], (*from).quals[i]);
        
        for(j=0;j<6;j++){
            (*target).base_counts[i][j]=(*from).base_counts[i][j];
            (*target).base_freqs[i][j]=(*from).base_freqs[i][j];
        }
        (*target).filtered_cov[i]=(*from).filtered_cov[i];
            
        (*target).ins_counts[i]=(*from).ins_counts[i];
        (*target).del_counts[i]=(*from).del_counts[i];
        
        (*target).ins_bases[i] = (char**) malloc( (*target).ins_counts[i] * sizeof(char*));
        (*target).del_bases[i] = (char**) malloc( (*target).del_counts[i] * sizeof(char*));
        for (j=0;j<(*target).ins_counts[i];i++){
            (*target).ins_bases[i][j] = (char*) malloc((strlen((*from).ins_bases[i][j])+1) * sizeof(char));
            strcpy( (*target).ins_bases[i][j], (*from).ins_bases[i][j]);
        } 
        for (j=0;j<(*target).del_counts[i];i++){
            (*target).del_bases[i][j] = (char*) malloc((strlen((*from).del_bases[i][j])+1) * sizeof(char));
            strcpy( (*target).del_bases[i][j], (*from).del_bases[i][j]);
        }
    }
    
    
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
    
    char c = bases[*base_ptr]; 
    if(quals[*qual_ptr] >= baseq_lim + 33 ){
        if(c=='.' || c==',' )      base_counts[REFBASE]++;
        else if(c=='A' || c=='a' ) base_counts[ABASE]++;
        else if(c=='C' || c=='c' ) base_counts[CBASE]++;
        else if(c=='G' || c=='g' ) base_counts[GBASE]++;
        else if(c=='T' || c=='t' ) base_counts[TBASE]++;
        (*filtered_cov)++;
    }
    
    //no Q filtering for deletions 
    if(c=='*' ) base_counts[DELBASE]++;
    
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
    (*base_ptr)+= offset-&bases[*base_ptr] + indel_len;
    return 0;
}

/* 
    parse an insertion from the bases and quals
*/
int handle_insertion(char* bases,int* ins_count,int* base_ptr){
    char* offset;
    (*ins_count)++;
    int indel_len=strtol(&bases[*base_ptr+1],&offset,10);
    (*base_ptr)+= offset-&bases[*base_ptr] + indel_len;
    return 0;
}


////////////////////////////////////////////////////////////////////////////
// Collect indels
////////////////////////////////////////////////////////////////////////////

/*
    collect indels in all samples
*/
int collect_indels_all_samples(struct Mpileup_line* my_pup_line){
    int i;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        collect_indels((*my_pup_line).bases[i],
                       &(*my_pup_line).ins_bases[i],
                       (*my_pup_line).ins_counts[i],
                       &(*my_pup_line).del_bases[i],
                       (*my_pup_line).del_counts[i]);
    }
    return 0;
}

/*
    collect the inserted, and deleted bases
*/
int collect_indels(char* bases, char*** ins_bases, int ins_count, char*** del_bases, int del_count){
    //free memory allocated before
    free_indel_bases(ins_bases,del_bases);
    //allocate new memory
    *ins_bases = (char**) malloc( (ins_count+1) * sizeof(char*));
    *del_bases = (char**) malloc( (del_count+1) * sizeof(char*));
    (*ins_bases)[ins_count] = (*del_bases)[del_count] = NULL;
   
    int i,del_c,ins_c; //pointers in data
    i = del_c = ins_c = 0;
    char* offset;
    while(bases[i]!=0){
        //beginning and end of the read signs
        if(bases[i] == '$' ) i++; //next
        else if(bases[i] == '^' ) i+=2; //jump next character (mapq too)
        //deletions
        else if(bases[i]=='-' ) {
            int indel_len=strtol(&bases[i+1],&offset,10);
            (*del_bases)[del_c] = (char*) malloc( (indel_len+1) * sizeof(char));
            memcpy((*del_bases)[del_c],offset,indel_len * sizeof(char));
            (*del_bases)[del_c][indel_len]=0;
            del_c++;
            i+= offset-&bases[i+1] + indel_len;
        }
        //insertions
        else if(bases[i]=='+' ){
            int indel_len=strtol(&bases[i+1],&offset,10);
            (*ins_bases)[ins_c] = (char*) malloc( (indel_len+1) * sizeof(char));
            memcpy((*ins_bases)[ins_c],offset,indel_len * sizeof(char));
            (*ins_bases)[ins_c][indel_len]=0;
            ins_c++;
            i+= offset-&bases[i+1] + indel_len;
        }
        //real base data
        else i++;
    }
    return 0;
}
   
/*
    free memory of indel bases
*/
int free_indel_bases(char*** ins_bases,char*** del_bases){
    int i=0;
    if (*ins_bases != NULL){
        while ( (*ins_bases)[i] != NULL ) { //last pointer should be null
            free((*ins_bases)[i]);
            i++;
        }
        free(*ins_bases);
    }
    i=0;
    if (*del_bases !=NULL){
        while ( (*del_bases)[i] != NULL ) { //last pointer should be null
            free((*del_bases)[i]);
            i++;
        }
        free(*del_bases);
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////
// Proximal gap filtering 
////////////////////////////////////////////////////////////////////////////
  
/*
    updates last indel gap position if there is one at the position
*/
int update_last_gap(struct Mpileup_line* my_pup_line, char** last_gap_chrom, int* last_gap_pos, int* is_gap){
    *is_gap=1;
    int i=0;
    for(i=0;i<(*my_pup_line).n_samples;i++){
        //last gap is either an insertion start deletion start, or a deleted base
        if ((*my_pup_line).ins_counts[i]!=0 || 
            (*my_pup_line).del_counts[i]!=0 || 
            (*my_pup_line).base_counts[i][DELBASE]!=0){
            //copy last chrom
            if ( *last_gap_chrom != NULL ) free(*last_gap_chrom);
            *last_gap_chrom = (char*) malloc( (strlen((*my_pup_line).chrom)+1) * sizeof(char));
            strcpy(*last_gap_chrom,(*my_pup_line).chrom);
            //pos
            *last_gap_pos=(*my_pup_line).pos;
            //set inidicator
            *is_gap=0;
            
            break;
        }
    }
    return 0;
}

/*
    if position is gap delete all potential mutations too close
*/
int proximal_gap_hindsight_filter(struct Mpileup_line* potential_mut_lines,int* mut_ptr,
                                  char* last_gap_chrom,int last_gap_pos,
                                 int proximal_gap_min_distance){
    int i; 
    for(i= *mut_ptr-1;i>=0;i--){
        //filter position if it is too close to last gap
        if (strcmp(last_gap_chrom,potential_mut_lines[i].chrom) == 0 && //same chrom
            last_gap_pos - potential_mut_lines[i].pos  < proximal_gap_min_distance ){ // proximal gap
            
            //delete this line
            //free resources
            free_mpileup_line(&(potential_mut_lines[i]));
            //decrement mut_pointer
            (*mut_ptr)--;
        }
        //if reached far enough we can break
        else break;
    }
    return 0;
        
}

/*
    prints and deletes mutations, which have survived the hindsight proximal gap filtering
*/
int flush_accepted_mutations(struct Mpileup_line* potential_mut_lines,
                             char* recent_chrom,int recent_pos,int* mut_ptr,
                             int proximal_gap_min_distance){
    int i,last_skipped; 
    //look for the first which is accepted by hindsight proximal gap filtering
    last_skipped=0;
    for(i = (*mut_ptr)-1;i>=0;i--){
        //skip position if it is too close to last gap
        if (strcmp(recent_chrom,potential_mut_lines[i].chrom) == 0 && //same chrom
            recent_pos - potential_mut_lines[i].pos  < proximal_gap_min_distance ){ // proximal gap possible
            continue;
        }
       else{
           last_skipped=i+1;
           break;
       }
    }
    //if all of them skipped nothing to flush
    if(last_skipped==0) return 0;
    
    //print the ones accepted
    for(i=0;i<last_skipped;i++){
        //print it and free resources
        print_mutation(&(potential_mut_lines[i]));
        free_mpileup_line(&(potential_mut_lines[i]));
    }
    
    //copy the last ones to the first places (preserve order)
    int new_mut_ptr=0;
    for(i=last_skipped;i < *mut_ptr;i++){
        //copy it into next empty place
        copy_mpileup_line(&(potential_mut_lines[new_mut_ptr]),&(potential_mut_lines[i]));
        //free resources
        free_mpileup_line(&(potential_mut_lines[i]));
        new_mut_ptr++;
    }
    *mut_ptr=new_mut_ptr;
    return 0;
} 
        

///////////////////////////////////////////////////////////////////////////
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
    print mutation
*/
int print_mutation(struct Mpileup_line* my_pup_line){
    printf("%d\t%s\t%d\t%c\t%c\n",(*my_pup_line).mut_sample_idx,(*my_pup_line).chrom,(*my_pup_line).pos,
                (*my_pup_line).ref_nuq,(*my_pup_line).mut_base);
    return 0;
}


/*
    calls mutation from frequencies
*/
int call_mutation(struct Mpileup_line* my_pup_line,double sample_mut_freq_limit,
                  double min_other_ref_freq_limit,int cov_limit,
                  char* last_gap_chrom, int last_gap_pos, int proximal_gap_min_distance){
    
    //filter position if it is too close to last gap
    if ( last_gap_chrom != NULL && //no gap yet
         strcmp(last_gap_chrom,(*my_pup_line).chrom) == 0 && //same chrom
         (*my_pup_line).pos - last_gap_pos < proximal_gap_min_distance ){ // proximal gap
        return 0;
    }
        
    
    double sample_mut_freq,min_other_ref_freq;
    int sample_idx;
    char mut_base;
    
    get_max_non_ref_freq(my_pup_line,&sample_mut_freq,&sample_idx,&mut_base);
    get_min_ref_freq(my_pup_line,&min_other_ref_freq,sample_idx);

    if (sample_mut_freq >= sample_mut_freq_limit && // mut freq larger than limit
        min_other_ref_freq > min_other_ref_freq_limit && //other sample ref_freq higher than limit 
        (*my_pup_line).filtered_cov[sample_idx] >= cov_limit ){ //coverage higher than limit
        
        (*my_pup_line).mut_base=mut_base;
        (*my_pup_line).mut_sample_idx=sample_idx;
        return 1;
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

