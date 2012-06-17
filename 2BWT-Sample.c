#include <stdio.h>
#include <stdlib.h>
#include "2BWT-Interface.h"

int main() {
    int i,j,k,c;

    //Variables for backward and forward search
    unsigned int l,r,rev_l,rev_r;
    //Variables for result
    unsigned int offset;
    int sequenceId;
    unsigned int saCount;
    
    //Variables for pattern
    char pattern[1024];
    strcpy(pattern,"AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA");
    int patternLength = strlen(pattern);

    // Load up the index with the below statement
    printf("Loading index ... "); 
    fflush(stdout);
    Idx2BWT * idx2BWT = BWTLoad2BWT("sampleRefSequence.fa.index",".sa");
    printf("DONE\n\n"); 
    
    // Convert the pattern into 2BWT recognised coding scheme
    BWTConvertPattern(idx2BWT,pattern,patternLength,pattern);
    
/* // The following performs a backward search of the pattern
 *     printf("Performing backward search of the pattern..\n");
 *     BWTSARangeInitial(idx2BWT,pattern[patternLength-1],&l,&r);
 *     for (i=patternLength-2;i>=0;i--) {
 *         BWTSARangeBackward(idx2BWT,pattern[i],&l,&r);
 *     }
 *     printf("SA Range being = %u %u (%u)\n\n",l,r,r-l+1);
 */

// The followin performs an exact extend search of the pattern
    int leav_len = patternLength;
    r = rev_r = idx2BWT->bwt->textLength;
    l = rev_l = 0;
    printf("Performing foreward exact extend of the pattern...\n");
    bwt_extend_exact(idx2BWT, pattern, patternLength, &leav_len, 0, &l, &r, &rev_l, &rev_r);

// The following output the first 5 position of the pattern
    printf("Reporting %d arbitrary occurrences..\n",j);
    for (i=l;i<=r;i++) {
        BWTRetrievePositionFromSAIndex(idx2BWT,i,&sequenceId,&offset);
        printf("Occurrence found in sequence #%d with offset %u\n",sequenceId,offset);
    }
    
    
    // Free up the 2BWT index
    printf("\nFree index ... "); 
    fflush(stdout);
    BWTFree2BWT(idx2BWT);
    printf("DONE\n"); 
    
    return 0;
}
