#include "isomut_lib.h"


////////////////////////////////////////////////////////////////////////////
// initialize
////////////////////////////////////////////////////////////////////////////
/*
    initializes all the pointers with NULL
        - always call this first after creating the struct,
            to avoid ERROR due to freeing/using memory garbage
*/
int init_mplp(struct mplp* my_mplp){
    int i;
    my_mplp->raw_line=NULL;
    my_mplp->chrom=NULL;
    my_mplp->pos=-42;
    my_mplp->n_samples=-42;
    my_mplp->ref_nuq='E';
    
    for(i=0;i<MAXSAMPLE;i++){
        my_mplp->raw_bases[i]=NULL;
        my_mplp->raw_quals[i]=NULL;
        my_mplp->ins_bases[i]=NULL;
        my_mplp->del_bases[i]=NULL;
    }
    
    
    strcpy(my_mplp->mut_type,"NOT\0");
    my_mplp->mut_base='E';
    my_mplp->mut_score=-42;
    for (i=0;i<MAX_INDEL_LEN;i++){my_mplp->mut_indel[i]='.';}
    my_mplp->mut_sample_idx=-42;
    my_mplp->mut_freq=-42;
    my_mplp->cleanliness=-42;
    
    return 0;
}

////////////////////////////////////////////////////////////////////////////
// free resources
////////////////////////////////////////////////////////////////////////////
/*
    frees all malloced objects in mplp struct
        - dont try to free uninitialized struct!
*/
int free_mplp(struct mplp* my_mplp){
    if( my_mplp->raw_line !=NULL) free(my_mplp->raw_line);
    if( my_mplp->chrom !=NULL) free(my_mplp->chrom);
    
    int i,j;
    for(i=0;i<my_mplp->n_samples;i++){
        if( my_mplp->raw_bases[i] !=NULL) free(my_mplp->raw_bases[i]);
        if( my_mplp->raw_quals[i] !=NULL) free(my_mplp->raw_quals[i]);
        
        if( my_mplp->ins_bases[i] !=NULL){
            for(j=0;j<my_mplp->counts[i][INS_START_IDX];j++){
                free(my_mplp->ins_bases[i][j]);
            }
            free(my_mplp->ins_bases[i]);
        }
        if( my_mplp->del_bases[i] !=NULL){
            for(j=0;j<my_mplp->counts[i][DEL_START_IDX];j++){
                free(my_mplp->del_bases[i][j]);
            }
            free(my_mplp->del_bases[i]);
        }
    }
    //initalize the pointers to null after freeing
    init_mplp(my_mplp);
    return 0;
}


////////////////////////////////////////////////////////////////////////////
// deep copy pileup struct
////////////////////////////////////////////////////////////////////////////

/*
    deep copy mpileup line
*/
int copy_mplp(struct mplp* target ,struct mplp* from) {
    //free and initialize resources int the target
    free_mplp(target);
    
    //copy raw line and chromosome
    target->raw_line = (char*) malloc((strlen(from->raw_line)+1) * sizeof(char));
    strcpy( target->raw_line, from->raw_line);
    target->chrom = (char*) malloc((strlen(from->chrom)+1) * sizeof(char));
    strcpy( target->chrom, from->chrom);
    //copy the position level data
    target->pos=from->pos;
    target->ref_nuq=from->ref_nuq;
    target->n_samples=from->n_samples;
    ///copy the mutation related data
    target->mut_base=from->mut_base;
    target->mut_sample_idx=from->mut_sample_idx;
    target->mut_score=from->mut_score;
    target->mut_freq=from->mut_freq;
    target->cleanliness=from->cleanliness;
    strcpy( target->mut_type, from->mut_type);
    strncpy( target->mut_indel, from->mut_indel,MAX_INDEL_LEN);
    //loop over sample level data
    int i,j;
    for(i=0;i<target->n_samples;i++){
        //copy basic data cov, bases, quals 
        target->raw_cov[i]=from->raw_cov[i];
        target->raw_bases[i] = (char*) malloc((strlen(from->raw_bases[i])+1) * sizeof(char));
        strcpy( target->raw_bases[i], from->raw_bases[i]);
        target->raw_quals[i] = (char*) malloc((strlen(from->raw_quals[i])+1) * sizeof(char));
        strcpy( target->raw_quals[i], from->raw_quals[i]);
        //copy filtered data, counts, coverage
        for(j=0;j<MAX_IDX;j++){
            target->counts[i][j]=from->counts[i][j];
            target->freqs[i][j]=from->freqs[i][j];
        }
        //copy all indel sequences
        target->ins_bases[i] = (char**) malloc(target->counts[i][INS_START_IDX] * sizeof(char*));
        target->del_bases[i] = (char**) malloc(target->counts[i][DEL_START_IDX] * sizeof(char*));
        for (j=0;j<target->counts[i][INS_START_IDX];j++){
            target->ins_bases[i][j] = (char*) malloc((strlen(from->ins_bases[i][j])+1) * sizeof(char));
            strcpy( target->ins_bases[i][j], from->ins_bases[i][j]);
        } 
        for (j=0;j<target->counts[i][DEL_START_IDX];j++){
            target->del_bases[i][j] = (char*) malloc((strlen(from->del_bases[i][j])+1) * sizeof(char));
            strcpy( target->del_bases[i][j], from->del_bases[i][j]);
        }
    }
    return 0;
}
    

////////////////////////////////////////////////////////////////////////////
// formatted printing functions
////////////////////////////////////////////////////////////////////////////


/*
    printf formatted representation of mpileup line
*/
int print_mplp(struct mplp* my_mplp){
    //print position level info
    printf("%s %d %c\n",my_mplp->chrom,my_mplp->pos,my_mplp->ref_nuq);
    
    //print counts and freqs and bases and quals
    int i;
    for(i=0;i<my_mplp->n_samples;i++){
        printf("cov %d,A %d,C %d,G %d,T %d,del %d,ins_start %d, del_start %d,read_start %d, read_end %d\n",
               my_mplp->counts[i][COV_IDX],
               my_mplp->counts[i][A_IDX],
               my_mplp->counts[i][C_IDX],
               my_mplp->counts[i][G_IDX],
               my_mplp->counts[i][T_IDX],
               my_mplp->counts[i][DEL_IDX],
               my_mplp->counts[i][INS_START_IDX],
               my_mplp->counts[i][DEL_START_IDX],
               my_mplp->counts[i][READ_START_IDX],
               my_mplp->counts[i][READ_END_IDX]);
        printf("A %.2f,C %.2f,G %.2f,T %.2f,del %.2f,ins_start %.2f, del_start %.2f,read_start %.2f, read_end %.2f\n",
               my_mplp->freqs[i][A_IDX],
               my_mplp->freqs[i][C_IDX],
               my_mplp->freqs[i][G_IDX],
               my_mplp->freqs[i][T_IDX],
               my_mplp->freqs[i][DEL_IDX],
               my_mplp->freqs[i][INS_START_IDX],
               my_mplp->freqs[i][DEL_START_IDX],
               my_mplp->freqs[i][READ_START_IDX],
               my_mplp->freqs[i][READ_END_IDX]);
        printf("%d %s %s\n",my_mplp->raw_cov[i],my_mplp->raw_bases[i],my_mplp->raw_quals[i]);
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////
// process input mpileup line from samtools mpileup command output
////////////////////////////////////////////////////////////////////////////

/*
      process input mpileup line from samtools mpileup command output
*/
int process_mplp_input_line(struct mplp* my_mplp,char* line, ssize_t line_size,int baseq_limit){
    
    //read mpileup 
    get_mplp(my_mplp,line,line_size);
    
    //count bases
    count_bases_all_samples(my_mplp,baseq_limit);
    
    //calculate freqs
    calculate_freqs_all_samples(my_mplp);
    
    //collect indels
    collect_indels_all_samples(my_mplp,baseq_limit);
    
    return 0;
}
    
    
    
////////////////////////////////////////////////////////////////////////////
// read pileup struct from mpileup line
////////////////////////////////////////////////////////////////////////////

/*
    gets mpileup line from the char* line
*/
int get_mplp(struct mplp* my_mplp,char* line, ssize_t line_size){   
    //store the raw line too
    free(my_mplp->raw_line);
    my_mplp->raw_line = (char*)malloc( (line_size+1) * sizeof(char));
    strcpy(my_mplp->raw_line,line);
    
    //temp buffer for reading those entries which will be formatted as not strings
    char* tmp_str=NULL;
    
    ssize_t i=0;
    while(i<line_size){
        //chrom
        get_next_entry(line,line_size,&i,&(my_mplp->chrom));
        //position
        get_next_entry(line,line_size,&i,&tmp_str);
        my_mplp->pos= (int) strtol(tmp_str,NULL,10);
        //ref nuq
        get_next_entry(line,line_size,&i,&tmp_str);
        my_mplp->ref_nuq=tmp_str[0];
        
        //read samples
        int temp_sample=0;
        while(i<line_size){
            //coverage
            get_next_entry(line,line_size,&i,&tmp_str);
            my_mplp->raw_cov[temp_sample]= (int) strtol(tmp_str,NULL,10);
            //bases
            get_next_entry(line,line_size,&i,&(my_mplp->raw_bases[temp_sample]));
            //quals
            get_next_entry(line,line_size,&i,&(my_mplp->raw_quals[temp_sample]));
            temp_sample++;
        }
        //store the number of samples
        my_mplp->n_samples=temp_sample;
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
    c0 = (int) *pointer;
    while(line[*pointer]!='\t' && *pointer<line_size)(*pointer)++;
    c = (int) *pointer-c0;
    (*pointer)++;
    
    if(*result != NULL) free(*result);
    *result = (char*)malloc( (c+1) * sizeof(char));
    memcpy(*result,line+c0,c * sizeof(char));
    (*result)[c] = 0;
    
    return c;
}



////////////////////////////////////////////////////////////////////////////
// Count bases in mplieup struct
////////////////////////////////////////////////////////////////////////////

/*
    counts bases,indels, and read start end signs in all samples
*/
int count_bases_all_samples(struct mplp* my_mplp,int baseq_lim){
    int i=0;
    for(i=0;i<my_mplp->n_samples;i++){
        count_bases(my_mplp->raw_bases[i],my_mplp->raw_quals[i],my_mplp->counts[i],
                    my_mplp->ref_nuq,baseq_lim);
    }
    return 0;
}

/*
    counts bases in one sample
*/
int count_bases(char* bases, char* quals,int* counts,char ref_base,int baseq_lim){
    //initialize counts to zero
    int i,j;
    for( i=0;i<MAX_IDX;i++) counts[i]=0;
    
    i = 0; //pointer in str for bases
    j = 0; //pointer in str for qualities
    counts[COV_IDX]=0;
    while(bases[i]!=0){
        //beginning and end of the read signs
        if(bases[i] == '^' ) {i+=2;counts[READ_START_IDX]++;} 
        else if(bases[i] == '$' ) {i++,counts[READ_END_IDX]++;} 
        //deletions
        else if(bases[i]=='-' ) handle_deletion(bases,&counts[DEL_START_IDX],&i,quals[j-1],baseq_lim);
        //insetions
        else if(bases[i]=='+' ) handle_insertion(bases,&counts[INS_START_IDX],&i,quals[j-1],baseq_lim); 
        //real base data
        else handle_base(bases,quals,counts,&counts[COV_IDX],&i,&j,baseq_lim);
    }
    //add refbase to corresponding base
    if(ref_base=='A') counts[A_IDX]+=counts[REF_IDX];
    else if(ref_base=='C') counts[C_IDX]+=counts[REF_IDX];
    else if(ref_base=='G') counts[G_IDX]+=counts[REF_IDX];
    else if(ref_base=='T') counts[T_IDX]+=counts[REF_IDX];
    return 0;
}


/* 
    parse a base from the bases and quals
*/
int handle_base(char* bases,char* quals,int* base_counts, int* filtered_cov,
                   int* base_ptr,int* qual_ptr,int baseq_lim){
    
    char c = bases[*base_ptr]; 
    if(quals[*qual_ptr] >= baseq_lim + 33 ){
        if(c=='.' || c==',' )      base_counts[REF_IDX]++;
        else if(c=='A' || c=='a' ) base_counts[A_IDX]++;
        else if(c=='C' || c=='c' ) base_counts[C_IDX]++;
        else if(c=='G' || c=='g' ) base_counts[G_IDX]++;
        else if(c=='T' || c=='t' ) base_counts[T_IDX]++;
        else if(c=='*' ) base_counts[DEL_IDX]++;
        (*filtered_cov)++;
    }
    (*qual_ptr)++;
    (*base_ptr)++;
    
    return 0;
}

/* 
    parse a deletion from the bases and quals
*/
int handle_deletion(char* bases,int* del_count,int* base_ptr,char qual,int baseq_lim){
    char* offset;
    int indel_len= (int) strtol(&bases[*base_ptr+1],&offset,10);
    (*base_ptr)+= offset-&bases[*base_ptr] + indel_len;
    if(qual >= baseq_lim + 33 )(*del_count)++;
    return 0;
}

/* 
    parse an insertion from the bases and quals
*/
int handle_insertion(char* bases,int* ins_count,int* base_ptr,char qual,int baseq_lim){
    char* offset;
    int indel_len= (int) strtol(&bases[*base_ptr+1],&offset,10);
    (*base_ptr)+= offset-&bases[*base_ptr] + indel_len;
    if(qual >= baseq_lim + 33 ) (*ins_count)++;
    return 0;
}


///////////////////////////////////////////////////////////////////////////
// Calculate base frequencies
////////////////////////////////////////////////////////////////////////////

/*
    calculate base freqs in all samples
*/
int calculate_freqs_all_samples(struct mplp* my_mplp){
    int i=0;
    for(i=0;i<my_mplp->n_samples;i++){
        calculate_freqs(my_mplp->freqs[i],my_mplp->counts[i]);
    }
    return 0;
}


/*
   calculate freqs in a sample
*/
int calculate_freqs(double* freqs,int* counts){
    int i;
    if (counts[COV_IDX]!=0){
        for(i=0;i<MAX_IDX;i++){
            freqs[i]=( (double) counts[i]) / counts[COV_IDX];
        }
    }
    else {
        for(i=0;i<MAX_IDX;i++){
            freqs[i]= ZERO_COV_FREQ;
        }
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////
// Collect indels
////////////////////////////////////////////////////////////////////////////

/*
    collect indels in all samples
*/
int collect_indels_all_samples(struct mplp* my_mplp,int baseq_lim){
    int i;
    for(i=0;i<my_mplp->n_samples;i++){
        collect_indels(my_mplp->raw_bases[i],
                       my_mplp->raw_quals[i],
                       &my_mplp->ins_bases[i],
                       my_mplp->counts[i][INS_START_IDX],
                       &my_mplp->del_bases[i],
                       my_mplp->counts[i][DEL_START_IDX],
                       baseq_lim);
    }
    return 0;
}

/*
    collect the inserted, and deleted bases
*/
int collect_indels(char* bases,char* quals, char*** ins_bases, int ins_count, 
                   char*** del_bases,int del_count, int baseq_lim){
    //allocate new memory
    if( *ins_bases!=NULL || *del_bases!=NULL){
        printf("ERROR: collect_indels() called on not NULL ins_bases, del_bases pointers,\n");
        printf("       maybe its called 2nd time, or mplp struct not freed, or not initialized,\n");
        printf("       possible memory leak, exiting");
        exit(1);
    }
    *ins_bases = (char**) malloc( (ins_count) * sizeof(char*));
    *del_bases = (char**) malloc( (del_count) * sizeof(char*));
   
    int i,j,del_c,ins_c; //pointers in data
    i = j = del_c = ins_c = 0;
    char* offset;
    while(bases[i]!=0){
        //beginning and end of the read signs
        if(bases[i] == '$' ) i++; //next
        else if(bases[i] == '^' ) i+=2; //jump next character (mapq too)
        //deletions
        else if(bases[i]=='-' ) {
            int indel_len= (int) strtol(&bases[i+1],&offset,10);
            i+= offset-&bases[i] + indel_len;
            if( quals[j-1] >= baseq_lim + 33 ){
                (*del_bases)[del_c] = (char*) malloc( (indel_len+1) * sizeof(char));
                memcpy((*del_bases)[del_c],offset,indel_len * sizeof(char));
                (*del_bases)[del_c][indel_len]=0;
                del_c++;
            }
        }
        //insertions
        else if(bases[i]=='+' ){
            int indel_len= (int) strtol(&bases[i+1],&offset,10);
            i+= offset-&bases[i] + indel_len;
            if( quals[j-1] >= baseq_lim + 33 ){
                (*ins_bases)[ins_c] = (char*) malloc( (indel_len+1) * sizeof(char));
                memcpy((*ins_bases)[ins_c],offset,indel_len * sizeof(char));
                (*ins_bases)[ins_c][indel_len]=0;
                ins_c++;
            }
        }
        //real base data
        else {i++;j++;}
    }
    return 0;
}
  

////////////////////////////////////////////////////////////////////////////
// Proximal gap filtering 
////////////////////////////////////////////////////////////////////////////
  
/*
    updates last indel gap position if there is one at the position
*/
int update_last_gap(struct mplp* my_mplp, char** last_gap_chrom, 
                    int* last_gap_pos_start, int* last_gap_pos_end, int* is_gap){
    *is_gap=1;
    int i,j;
    for(i=0;i<my_mplp->n_samples;i++){
        //last gap is either an insertion start deletion start or 0 cov
        // or more than 50% read start, read end
        if (my_mplp->counts[i][INS_START_IDX]!=0 || 
            my_mplp->counts[i][DEL_START_IDX]!=0 ||
            my_mplp->freqs[i][READ_START_IDX] > 0.5 ||
            my_mplp->freqs[i][READ_END_IDX] > 0.5 ){
            
            //copy last chrom 
            if ( *last_gap_chrom != NULL ) free(*last_gap_chrom);
            *last_gap_chrom = (char*) malloc( (strlen(my_mplp->chrom)+1) * sizeof(char));
            strcpy(*last_gap_chrom,my_mplp->chrom);
            //update last gap pos start
            *last_gap_pos_start = my_mplp->pos;
            //update last gap pos end for 0 coverage, read starts, ends
            if ( *last_gap_pos_start > *last_gap_pos_end ){
                 *last_gap_pos_end = *last_gap_pos_start;
             }
            //update last gap pos end for indels 
             if ( *last_gap_pos_start >= *last_gap_pos_end  && my_mplp->counts[i][INS_START_IDX]!=0 ){
                 *last_gap_pos_end = *last_gap_pos_start + 1;
             }
            //gap pos end is hihger for deletions
            for (j=0;j < my_mplp->counts[i][DEL_START_IDX];j++){
                if(*last_gap_pos_end < my_mplp->pos + 1 + (int) strlen(my_mplp->del_bases[i][j])){
                    *last_gap_pos_end = my_mplp->pos + 1 + (int) strlen(my_mplp->del_bases[i][j]);
                }
            }
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
int proximal_gap_hindsight_filter(struct mplp* saved_mutations,int* mut_ptr,
                                  char* last_gap_chrom,int last_gap_pos_start,
                                 int proximal_gap_min_dist_SNV,int proximal_gap_min_dist_indel){
    
    //if position is a called indel the filter has been run
    if( *mut_ptr > 0 &&  strcmp(last_gap_chrom,saved_mutations[(*mut_ptr) -1].chrom) == 0 &&
            last_gap_pos_start == saved_mutations[(*mut_ptr) -1].pos ){
        return 0;
    }
    int i,j; 
    for(i= *mut_ptr-1;i>=0;i--){
        //delet SNV too close to gap
        if (strcmp(last_gap_chrom,saved_mutations[i].chrom) != 0) return 0;
        else if( strcmp(saved_mutations[i].mut_type,"SNV")==0 &&
                    last_gap_pos_start - saved_mutations[i].pos  < proximal_gap_min_dist_SNV ){
            //delete this line
            free_mplp(&saved_mutations[i]);
            (*mut_ptr)--;
        }
        //delete indels too close to gap and quit, because filter was run then
        else if( ( strcmp(saved_mutations[i].mut_type,"INS")==0 && //ins
                   last_gap_pos_start - saved_mutations[i].pos  < proximal_gap_min_dist_indel)
                   || ( strcmp(saved_mutations[i].mut_type,"DEL")==0  &&  //del
                   last_gap_pos_start - saved_mutations[i].pos -  
                   ((int)strlen(saved_mutations[i].mut_indel))  < proximal_gap_min_dist_indel) ){
            //delete this line
            free_mplp(&saved_mutations[i]);
            (*mut_ptr)--;
            //copy all later SNVs one position ahead
            for(j=i;j < *mut_ptr;j++){
                copy_mplp(&saved_mutations[j],&saved_mutations[j+1]);
                free_mplp(&saved_mutations[j+1]);
            }
            return 0;
        }
    }
    return 0;
        
}

/*
    prints and deletes mutations, which have survived the hindsight proximal gap filtering
*/
int flush_accepted_mutations(struct mplp* saved_mutations,
                             char* recent_chrom,int recent_pos,int* mut_ptr,
                             int proximal_gap_min_dist_SNV,int proximal_gap_min_dist_indel){
    int i,j; 
    for(i = 0; i<*mut_ptr;i++){
        //check if mut is accepted already
        if ( ( strcmp(recent_chrom,saved_mutations[i].chrom) != 0 ) || //other chrom
             ( strcmp(saved_mutations[i].mut_type,"SNV")==0 && //SNV
             recent_pos - saved_mutations[i].pos  > proximal_gap_min_dist_SNV ) ||
             ((strcmp(saved_mutations[i].mut_type,"INS")==0 || //indel
               strcmp(saved_mutations[i].mut_type,"DEL")==0 ) && 
             recent_pos - saved_mutations[i].pos  > proximal_gap_min_dist_indel )){
            //print and delete it
            print_mutation(&saved_mutations[i]);
            free_mplp(&saved_mutations[i]);
            (*mut_ptr)--;
            //step everyone else ahead
            for(j=i;j < *mut_ptr;j++){
                copy_mplp(&saved_mutations[j],&saved_mutations[j+1]);
                free_mplp(&saved_mutations[j+1]);
            } 
            i--;
        }
    }
    return 0;
} 
        


////////////////////////////////////////////////////////////////////////////
// Call Mutations 
////////////////////////////////////////////////////////////////////////////

/*
    print mutation
*/
int print_mutation(struct mplp* my_mplp){
    if ( strcmp(my_mplp->mut_type,"SNV") ==0 ){
        printf("%d\t%s\t%d\t%s\t%.2f\t%c\t%c\t%d\t%.3f\t%.3f\n",
               my_mplp->mut_sample_idx,
               my_mplp->chrom,
               my_mplp->pos,
               my_mplp->mut_type,
               my_mplp->mut_score,
               my_mplp->ref_nuq,
               my_mplp->mut_base,
               my_mplp->counts[my_mplp->mut_sample_idx][COV_IDX],
               my_mplp->mut_freq,
               my_mplp->cleanliness);
    }
    if ( strcmp(my_mplp->mut_type,"INS") ==0   ){
        printf("%d\t%s\t%d\t%s\t%.2f\t-\t%s\t%d\t%.3f\t%.3f\n",
               my_mplp->mut_sample_idx,
               my_mplp->chrom,
               my_mplp->pos,
               my_mplp->mut_type,
               my_mplp->mut_score,
               my_mplp->mut_indel,
               my_mplp->counts[my_mplp->mut_sample_idx][COV_IDX],
               my_mplp->mut_freq,
               my_mplp->cleanliness);

    }
    if ( strcmp(my_mplp->mut_type,"DEL") ==0   ){
        printf("%d\t%s\t%d\t%s\t%.2f\t%s\t-\t%d\t%.3f\t%.3f\n",
               my_mplp->mut_sample_idx,
               my_mplp->chrom,
               my_mplp->pos,
               my_mplp->mut_type,
               my_mplp->mut_score,
               my_mplp->mut_indel,
               my_mplp->counts[my_mplp->mut_sample_idx][COV_IDX],
               my_mplp->mut_freq,
               my_mplp->cleanliness);

    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////
// Call SNVs
////////////////////////////////////////////////////////////////////////////


/*
    calls mutation from frequencies
*/
int call_snv(struct mplp* saved_mutations, int* mut_ptr, struct mplp* my_mplp,
             double sample_mut_freq_limit,double min_other_ref_freq_limit,int cov_limit,
             char* last_gap_chrom, int last_gap_pos_end, int proximal_gap_min_distance){
    
    //filter position if it is too close to last gap
    if ( last_gap_chrom != NULL && //no gap yet
         strcmp(last_gap_chrom,my_mplp->chrom) == 0 && //same chrom
         my_mplp->pos - last_gap_pos_end < proximal_gap_min_distance ){ // proximal gap
        return 0;
    }
        
    double sample_mut_freq,min_other_ref_freq;
    int sample_idx,other_idx;
    char mut_base = 'E';
    
    get_max_non_ref_freq(my_mplp,&sample_mut_freq,&sample_idx,&mut_base);
    get_min_ref_freq(my_mplp,&min_other_ref_freq,sample_idx,&other_idx);

    if (sample_mut_freq >= sample_mut_freq_limit && // mut freq larger than limit
        min_other_ref_freq > min_other_ref_freq_limit && //other sample ref_freq higher than limit 
        my_mplp->counts[sample_idx][COV_IDX] >= cov_limit ){ //coverage higher than limit
        
        my_mplp->mut_base=mut_base;
        my_mplp->mut_sample_idx=sample_idx;
        my_mplp->mut_freq=sample_mut_freq;
        my_mplp->cleanliness=min_other_ref_freq;
        strncpy(my_mplp->mut_type, "SNV\0",4);
       
        my_mplp->mut_score = -log10(fisher22((uint32_t) ((1-sample_mut_freq) * my_mplp->counts[sample_idx][COV_IDX]),
                               (uint32_t) (sample_mut_freq * my_mplp->counts[sample_idx][COV_IDX]),
                               (uint32_t) (min_other_ref_freq * my_mplp->counts[other_idx][COV_IDX]),
                               (uint32_t) ((1-min_other_ref_freq) * my_mplp->counts[other_idx][COV_IDX]),1));
        
        //save potential mutation
        copy_mplp(&saved_mutations[*mut_ptr],my_mplp);
        (*mut_ptr)++;
    }
    return 0;        
}


/*
    gets the highest not reference mut freq
*/
int get_max_non_ref_freq(struct mplp* my_mplp,double* val, int* idx, char* mut_base){

    //arrays only for encoding/decoding base indices to bases
    int base_2_idx[256];
    base_2_idx[ 'A' ] = A_IDX;base_2_idx[ 'C' ] = C_IDX;
    base_2_idx[ 'G' ] = G_IDX;base_2_idx[ 'T' ] = T_IDX;
    char idx_2_base[MAX_IDX];
    idx_2_base[ A_IDX ] = 'A';idx_2_base[ C_IDX ] = 'C';
    idx_2_base[ G_IDX ] = 'G';idx_2_base[ T_IDX ] = 'T';
  
    //get index for ref base
    int ref_idx=base_2_idx[(int)my_mplp->ref_nuq];
   
    //initialize values to negative numbers
    *val = -42;
    *idx  = -42;
    
    int i,j;
    //loop over samples
    for(i=0;i<my_mplp->n_samples;i++){
        //loop over bases
        for(j=A_IDX;j<T_IDX+1;j++){
            if ( j!=ref_idx && //base is not the reference base
                my_mplp->freqs[i][j] > *val && //larger than largest yet
                my_mplp->freqs[i][j] != ZERO_COV_FREQ ){ //not a 0 cov sample
                //save value of max, sample idx, and the mut base as chr
                *val=my_mplp->freqs[i][j];
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
int get_min_ref_freq(struct mplp* my_mplp,double* val, int idx_2skip, int* other_idx ){
    //initialize value to large number
    *val = 42;
    *other_idx=0;
    
    int i;
    //loop over samples
    for(i=0;i<my_mplp->n_samples;i++){
        if ( i != idx_2skip && // skip the mutated sample
           my_mplp->freqs[i][REF_IDX] < *val && //smaller than smallest yet
           my_mplp->freqs[i][REF_IDX] != ZERO_COV_FREQ ){ //not a 0 cov sample
                //save value of min
            *val=my_mplp->freqs[i][REF_IDX];
            *other_idx=i;
        }
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////
// Call indels 
////////////////////////////////////////////////////////////////////////////


/*
    calls indels
*/
int call_indel(struct mplp* saved_mutations, int* mut_ptr, struct mplp* my_mplp,
               double sample_mut_freq_limit,double min_other_ref_freq_limit,int cov_limit,
               char* last_gap_chrom, int last_gap_pos_start,int last_gap_pos_end,
               int prox_gap_min_dist_SNV,int prox_gap_min_dist_indel){
    
    //filter position if it is too close to last gap
    if ( last_gap_chrom != NULL && //no gap yet
         strcmp(last_gap_chrom,my_mplp->chrom) == 0  && //same chrom
         my_mplp->pos - last_gap_pos_end < prox_gap_min_dist_indel  && // proximal gap
         my_mplp->pos != last_gap_pos_start ){ //  dont filter for indel because of himself!
        return 0;
    }
    
    
    double sample_indel_freq,min_other_noindel_freq;
    int sample_idx,other_idx;
    char mut_indel[MAX_INDEL_LEN];
    char mut_type[4];
    
    get_max_indel_freq(my_mplp,&sample_indel_freq,&sample_idx,mut_indel,mut_type);
    get_min_other_noindel_freq(my_mplp,&min_other_noindel_freq,sample_idx,&other_idx);
    
    //todo??
    //filter for non unique indels? 
    
    if (sample_indel_freq >= sample_mut_freq_limit && // indel freq larger than limit
        min_other_noindel_freq > min_other_ref_freq_limit && //  cleanness 
        my_mplp->counts[sample_idx][COV_IDX] >= cov_limit ){ //coverage higher than limit
        
        strncpy(my_mplp->mut_indel,mut_indel,MAX_INDEL_LEN);
        my_mplp->mut_sample_idx=sample_idx;
        my_mplp->mut_freq=sample_indel_freq;
        my_mplp->cleanliness=min_other_noindel_freq;
        strncpy(my_mplp->mut_type,mut_type,4);
       
        my_mplp->mut_score = -log10(fisher22((uint32_t) ((1-sample_indel_freq) * my_mplp->counts[sample_idx][COV_IDX]),
                               (uint32_t) (sample_indel_freq * my_mplp->counts[sample_idx][COV_IDX]),
                               (uint32_t) (min_other_noindel_freq * my_mplp->counts[other_idx][COV_IDX]),
                               (uint32_t) ((1-min_other_noindel_freq) *my_mplp->counts[other_idx][COV_IDX]),1));
        
        //proximal hindsight filtering for SNVs before
        proximal_gap_hindsight_filter(saved_mutations,mut_ptr,my_mplp->chrom,
                                      my_mplp->pos,prox_gap_min_dist_SNV,prox_gap_min_dist_indel);
        
        //save potential mutation
        copy_mplp(&saved_mutations[*mut_ptr],my_mplp);
        (*mut_ptr)++;
    }
    return 0;        
}



/*
    gets the highest indel freq 
*/
int get_max_indel_freq(struct mplp* my_mplp,double* val, int* idx, char* mut_indel,char* mut_type){
    //initialize values to negative numbers
    *val = 0 ;
    *idx  = 0 ; 
    
    int i;
    //loop over samples
    for(i=0;i<my_mplp->n_samples;i++){
        if( my_mplp->freqs[i][INS_START_IDX] > (*val) && //larger than largest yet
            my_mplp->freqs[i][INS_START_IDX] != ZERO_COV_FREQ ){ //not a 0 cov sample
           //save value of max, sample idx, and the mut base as chr
            *val=my_mplp->freqs[i][INS_START_IDX];
            *idx=i;
            strncpy(mut_indel,my_mplp->ins_bases[i][0],strlen(my_mplp->ins_bases[i][0])+1);
            strncpy(mut_type,"INS\0",4);
        }
        if( my_mplp->freqs[i][DEL_START_IDX] > (*val) && //larger than largest yet
            my_mplp->freqs[i][DEL_START_IDX] != ZERO_COV_FREQ ){ //not a 0 cov sample
           //save value of max, sample idx, and the mut base as chr
            *val=my_mplp->freqs[i][DEL_START_IDX];
            *idx=i;
            strncpy(mut_indel,my_mplp->del_bases[i][0],strlen(my_mplp->del_bases[i][0])+1);
            strncpy(mut_type,"DEL\0",4);
        }
    }
    return 0;
}


/*
    gets the highest indel freq, except for 1 sample
*/
int get_min_other_noindel_freq(struct mplp* my_mplp,double* val, int idx_2skip, int* other_idx ){
    //initialize value to  large number
    *val =  42;
    *other_idx=0;
    
    int i;
    //loop over samples
    for(i=0;i<my_mplp->n_samples;i++){
        double noindel_freq= my_mplp->freqs[i][REF_IDX] - my_mplp->freqs[i][INS_START_IDX] - 
                               my_mplp->freqs[i][DEL_START_IDX];
        if ( i != idx_2skip && // skip the mutated sample
               noindel_freq < *val  && // smaller 
               my_mplp->freqs[i][INS_START_IDX] != ZERO_COV_FREQ ){ //not a 0 cov sample
            //save value of min 
            *val=noindel_freq;
            *other_idx=i;
        }
    }
    return 0;
}
