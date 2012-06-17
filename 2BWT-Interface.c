/*

   Date   : 19th June 2011
   Author : Edward MK Wu
   Change : Packaging 2BWT library as a separate product.
            Enhancing the 2BWT library with interface functions
            for basic BWT search.

*/

#include "2BWT-Interface.h"

Idx2BWT * BWTLoad2BWT(const char * indexFilePrefix, const char * saFileNameExtension) {

    Idx2BWT * idx2BWT = malloc(sizeof(Idx2BWT));
    BWT * bwt;
    BWT * rev_bwt;
    HSP * hsp;
    
    char bwtFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char bwtOccFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char saFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char rev_bwtFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char rev_bwtOccFilename[MAX_INDEX_FILENAME_LENGTH];
    char packedDnaFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char annotationFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char ambiguityFilename[MAX_INDEX_FILENAME_LENGTH]; 
    char translateFilename[MAX_INDEX_FILENAME_LENGTH]; 
    
    strcpy(bwtFilename,indexFilePrefix);
    strcpy(bwtOccFilename,indexFilePrefix);
    strcpy(saFilename,indexFilePrefix);
    strcpy(rev_bwtFilename,indexFilePrefix);
    strcpy(rev_bwtOccFilename,indexFilePrefix);
    strcpy(packedDnaFilename,indexFilePrefix);
    strcpy(annotationFilename,indexFilePrefix);
    strcpy(ambiguityFilename,indexFilePrefix);
    strcpy(translateFilename,indexFilePrefix);
    
    strcat(bwtFilename,".bwt");
    strcat(bwtOccFilename,".fmv");
    strcat(saFilename,saFileNameExtension);
    strcat(rev_bwtFilename,".rev.bwt");
    strcat(rev_bwtOccFilename,".rev.fmv");
    strcat(packedDnaFilename,".pac");
    strcat(annotationFilename,".ann");
    strcat(ambiguityFilename,".amb");
    strcat(translateFilename,".tra");
    
    MMMasterInitialize(3, 0, FALSE, NULL);
    MMPool * mmPool = MMPoolCreate(2097152);
    
    bwt = BWTLoad(mmPool, bwtFilename, bwtOccFilename, saFilename, NULL, NULL, NULL);
    rev_bwt = BWTLoad(mmPool, rev_bwtFilename, rev_bwtOccFilename, NULL, NULL, NULL, NULL);
    hsp = HSPLoad(mmPool, packedDnaFilename, annotationFilename, ambiguityFilename,translateFilename, 1);
    
    HSPFillCharMap(idx2BWT->charMap);
    HSPFillComplementMap(idx2BWT->complementMap);
     
    idx2BWT->bwt = bwt;
    idx2BWT->rev_bwt = rev_bwt;
    idx2BWT->hsp = hsp;
    idx2BWT->mmPool = mmPool;
    
    return idx2BWT;
}

void BWTFree2BWT(Idx2BWT * idx2BWT) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    HSP * hsp = idx2BWT->hsp;
    MMPool * mmPool = idx2BWT->mmPool;
    
    HSPFree(mmPool, hsp, 1);
    BWTFree(mmPool, bwt);
    BWTFree(mmPool, rev_bwt);
    MMPoolFree(mmPool);
    
    free(idx2BWT);
}

void BWTConvertPattern(Idx2BWT * idx2BWT, const char * patternSource, int patternLength, unsigned char * patternDestination) {

    int i;
    for (i=0;i<patternLength;i++) {
        patternDestination[i] = idx2BWT->charMap[patternSource[i]];
    }
    patternDestination[i]='\0';
}


void BWTSARangeInitial(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    //TODO I guess saIndex is start from 1, and cumulativeFreq is from 0
    (*saIndexLeft) = bwt->cumulativeFreq[c]+1;
    (*saIndexRight) = bwt->cumulativeFreq[c+1];

}


//给出saIndex, backward 算出 saIndex
void BWTSARangeBackward(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned int l = (*saIndexLeft);
    unsigned int r = (*saIndexRight);
    (*saIndexLeft)  = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c) + 1;
    (*saIndexRight) = bwt->cumulativeFreq[c] + BWTOccValue(bwt, r + 1, c);

}

//给出saIndex, foreward 算出 saIndex
void BWTSARangeForeward(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight) {

    BWT * rev_bwt = idx2BWT->rev_bwt;
    BWT *bwt = idx2BWT->bwt;
    
    unsigned int l = (*saIndexLeft);
    unsigned int r = (*saIndexRight);
    (*saIndexLeft)  = bwt->cumulativeFreq[c] + BWTOccValue(rev_bwt, l, c) + 1;
    (*saIndexRight) = bwt->cumulativeFreq[c] + BWTOccValue(rev_bwt, r + 1, c);

}

//给出saIndex,rev_saIndex backward算出saIndex 和 rev_saIndex
void BWTSARangeBackward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight,
                        unsigned int *rev_saIndexLeft, unsigned int *rev_saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    unsigned int l = (*saIndexLeft);
    unsigned int r = (*saIndexRight);
    unsigned int rev_l = (*rev_saIndexLeft);
    unsigned int rev_r = (*rev_saIndexRight);
    int k;
    
    BWTAllOccValue(bwt,l,oL);
    BWTAllOccValue(bwt,r + 1,oR); //TODO why r need plus 1
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }

    l = bwt->cumulativeFreq[c] + oL[c] + 1;
    r = bwt->cumulativeFreq[c] + oR[c];
    rev_r = rev_r - oCount[c];
    rev_l = rev_r - (r-l);

    (*saIndexLeft) = l;
    (*saIndexRight) = r;
    (*rev_saIndexLeft) = rev_l;
    (*rev_saIndexRight) = rev_r;
    
}



//根据rev_saIndex saIndex, forward算出rev_saIndex and saIndex
void BWTSARangeForward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight,
                        unsigned int *rev_saIndexLeft, unsigned int *rev_saIndexRight) {

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    unsigned int l = (*saIndexLeft);
    unsigned int r = (*saIndexRight);
    unsigned int rev_l = (*rev_saIndexLeft);
    unsigned int rev_r = (*rev_saIndexRight);
    int k;
    
    BWTAllOccValue(rev_bwt,rev_l,oL);
    BWTAllOccValue(rev_bwt,rev_r + 1,oR);
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }
    
    rev_l = bwt->cumulativeFreq[c] + oL[c] + 1;
    rev_r = bwt->cumulativeFreq[c] + oR[c];
    r = r - oCount[c];
    l = r - (rev_r-rev_l);
    
    (*saIndexLeft) = l;
    (*saIndexRight) = r;
    (*rev_saIndexLeft) = rev_l;
    (*rev_saIndexRight) = rev_r;
}


//根据给出的saIndex, backward算出所有字符的saIndex
void BWTAllSARangesBackward(Idx2BWT * idx2BWT, const unsigned int saIndexLeft,
        const unsigned int saIndexRight,unsigned int *resultSaIndexesLeft,
        unsigned int *resultSaIndexesRight)
{

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    unsigned int l = saIndexLeft;
    unsigned int r = saIndexRight;
    int k;
    
    BWTAllOccValue(bwt,l,oL);
    BWTAllOccValue(bwt,r + 1,oR);
    
    for (k=0;k<ALPHABET_SIZE;k++) {
        resultSaIndexesLeft[k]  = bwt->cumulativeFreq[k] + oL[k] + 1;
        resultSaIndexesRight[k] = bwt->cumulativeFreq[k] + oR[k];
    }

}

//根据给出的saIndex rev_saIndex, backward算出所有char的saIndex 和 rev_saIndex
void BWTAllSARangesBackward_Bidirection(Idx2BWT * idx2BWT, const unsigned int saIndexLeft,
        const unsigned int saIndexRight, const unsigned int rev_saIndexLeft,
        const unsigned int rev_saIndexRight, unsigned int *resultSaIndexLeft,
        unsigned int *resultSaIndexRight, unsigned int *rev_resultSaIndexLeft,
        unsigned int *rev_resultSaIndexRight)
{

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    unsigned int l = saIndexLeft;
    unsigned int r = saIndexRight;
    unsigned int rev_l = rev_saIndexLeft;
    unsigned int rev_r = rev_saIndexRight;
    int k;
    
    BWTAllOccValue(bwt,l,oL);
    BWTAllOccValue(bwt,r + 1,oR);
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }
    
    for (k=0;k<ALPHABET_SIZE;k++) {

        resultSaIndexLeft[k] = bwt->cumulativeFreq[k] + oL[k] + 1;
        resultSaIndexRight[k] = bwt->cumulativeFreq[k] + oR[k];
        rev_resultSaIndexRight[k] = rev_r - oCount[k];
        rev_resultSaIndexLeft[k] = rev_resultSaIndexRight[k] - 
                                    (resultSaIndexRight[k]-resultSaIndexLeft[k]);

    }

}

//根据rev_saIndex saIndex, forward算出所有char的saIndex and rev_saIndex
void BWTAllSARangesForward_Bidirection(Idx2BWT * idx2BWT, const unsigned int saIndexLeft,
        const unsigned int saIndexRight,const unsigned int rev_saIndexLeft,
        const unsigned int rev_saIndexRight,unsigned int *resultSaIndexLeft,
        unsigned int *resultSaIndexRight,unsigned int *rev_resultSaIndexLeft,
        unsigned int *rev_resultSaIndexRight)
{

    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    unsigned int l = saIndexLeft;
    unsigned int r = saIndexRight;
    unsigned int rev_l = rev_saIndexLeft;
    unsigned int rev_r = rev_saIndexRight;
    int k;
    
    BWTAllOccValue(rev_bwt,rev_l,oL);
    BWTAllOccValue(rev_bwt,rev_r + 1,oR);
    oCount[ALPHABET_SIZE-1]=0;
    for (k=ALPHABET_SIZE-2;k>=0;k--) {
        oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
    }
    for (k=0;k<ALPHABET_SIZE;k++) {
        rev_resultSaIndexLeft[k] = bwt->cumulativeFreq[k] + oL[k] + 1;
        rev_resultSaIndexRight[k] = bwt->cumulativeFreq[k] + oR[k];
        resultSaIndexRight[k] = r - oCount[k];
        resultSaIndexLeft[k] = resultSaIndexRight[k] - (rev_resultSaIndexRight[k]-rev_resultSaIndexLeft[k]);
    }
}

//根据saIndex还原出seq在ref中的位置
void BWTRetrievePositionFromSAIndex(Idx2BWT * idx2BWT, unsigned int saIndex,
        unsigned int * sequenceId, unsigned int * offset)
{
    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    HSP * hsp = idx2BWT->hsp;
	unsigned short * ambiguityMap = hsp->ambiguityMap;
	Translate * translate = hsp->translate;
    
    unsigned int ambPosition = BWTSaValue(bwt,saIndex);
    unsigned int approxIndex = ambPosition>>GRID_SAMPLING_FACTOR_2_POWER;
    unsigned int approxValue = ambiguityMap[approxIndex];
    while (translate[approxValue].startPos>ambPosition) {
        approxValue--;
    }
    ambPosition-=translate[approxValue].correction;
    
    (*sequenceId) = translate[approxValue].chrID;
    (*offset) = ambPosition;
}

//exact match seq to ref
int bwt_match_exact(Idx2BWT *bi_bwt, ubyte_t *seq, int len, bwtint_t *sa_begin,
        bwtint_t *sa_end, bwtint_t *rev_sa_begin, bwtint_t *rev_sa_end)
{
    bwtint_t l, k, rev_l, rev_k;
    int i;
    k = *sa_begin;
    l = *sa_end;
    rev_k = *rev_sa_begin;
    rev_l = *rev_sa_end;
    for (i = len - 1; i >= 0; i--) {
        if(seq[i]>3)
            return 0;
        BWTSARangeBackward_Bidirection(bi_bwt, seq[i], &k, &l, &rev_k, &rev_l);
        if(k > l)
            break;  //no match
    }
    if (k > l)
        return 0;
    if(*sa_begin) *sa_begin = k;
    if(*sa_end) *sa_end = l;
    if(*rev_sa_begin) *rev_sa_begin = rev_k;
    if(*rev_sa_end) *rev_sa_end = rev_l;
    return l - k + 1;
}

// exact extend
// leav_len: leaving length to extend
// type: 1 - backward, 0 - forward
int bwt_extend_exact(const Idx2BWT *bi_bwt, ubyte_t *seq, int len, int *leav_len, int type,
        bwtint_t *sa_begin, bwtint_t *sa_end, bwtint_t *rev_sa_begin, bwtint_t *rev_sa_end)
{
    bwtint_t k, l, rev_k, rev_l;
    int i;
    k = *sa_begin;
    l = *sa_end;
    rev_k = *rev_sa_begin;
    rev_l = *rev_sa_end;
    if(type == 1){
        for(i = *leav_len; i > 0; --i){
            if(seq[i - 1] > 3)
                break;
            BWTSARangeBackward_Bidirection(bi_bwt, seq[i], &k, &l, &rev_k, &rev_l);
            if(k > l)
                break;
            else{
                *sa_begin = k;
                *sa_end = l;
                *rev_sa_begin = rev_k;
                *rev_sa_end = rev_l;
            }
        }
    } else {
        for(i = (*leav_len); i > 0; --i){
            if(seq[len-i] > 3)
                break; 
            BWTSARangeForward_Bidirection(bi_bwt, seq[len-i], &k, &l, &rev_k, &rev_l);
            if(k > l)
                break;
            else{
                *sa_begin = k;
                *sa_end = l;
                *rev_sa_begin = rev_k;
                *rev_sa_end = rev_l;
            }
        }
    }
    *leav_len = i;
    return *sa_end - (*sa_begin) + 1;
}
