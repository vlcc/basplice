/*

   Date   : 19th June 2011
   Author : Edward MK Wu
   Change : Packaging 2BWT library as a separate product.
            Enhancing the 2BWT library with interface functions
            for basic BWT search.

*/

#ifndef __2BWT_INTERFACE_H__
#define __2BWT_INTERFACE_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "BWT.h"
#include "BWTLoad.h"
#include "HSP.h"
#include "MemManager.h"

#define MAX_INDEX_FILENAME_LENGTH 1024
#define POOLSIZE                  2097152

typedef uint32_t bwtint_t;
typedef unsigned char ubyte_t;

typedef struct _Idx2BWT{
    MMPool * mmPool;
    BWT * bwt;
    BWT * rev_bwt;
    HSP * hsp;
    unsigned char charMap[256];
    unsigned char complementMap[256];
} Idx2BWT ;

Idx2BWT *BWTLoad2BWT(const char * indexFilePrefix, const char * saFileNameExtension);

void BWTFree2BWT(Idx2BWT * idx2BWT);

void BWTConvertPattern(Idx2BWT * idx2BWT, const char * patternSource, int patternLength, unsigned char * patternDestination);
void BWTSARangeInitial(Idx2BWT * idx2BWT, const unsigned char c, unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight);

void BWTSARangeBackward(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight);

void BWTSARangeForeward(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight);

void BWTSARangeBackward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight,
                        unsigned int *rev_saIndexLeft, unsigned int *rev_saIndexRight);

void BWTSARangeForward_Bidirection(Idx2BWT * idx2BWT, const unsigned char c, 
                        unsigned int *saIndexLeft, unsigned int *saIndexRight,
                        unsigned int *rev_saIndexLeft, unsigned int *rev_saIndexRight);
                        
void BWTAllSARangesBackward(Idx2BWT * idx2BWT, 
                        const unsigned int saIndexLeft, const unsigned int saIndexRight,
                        unsigned int *resultSaIndexesLeft, unsigned int *resultSaIndexesRight);

void BWTAllSARangesBackward_Bidirection(Idx2BWT * idx2BWT, 
                        const unsigned int saIndexLeft, const unsigned int saIndexRight,
                        const unsigned int rev_saIndexLeft, const unsigned int rev_saIndexRight,
                        unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight,
                        unsigned int *rev_resultSaIndexLeft, unsigned int *rev_resultSaIndexRight);

void BWTAllSARangesForward_Bidirection(Idx2BWT * idx2BWT, 
                        const unsigned int saIndexLeft, const unsigned int saIndexRight,
                        const unsigned int rev_saIndexLeft, const unsigned int rev_saIndexRight,
                        unsigned int *resultSaIndexLeft, unsigned int *resultSaIndexRight,
                        unsigned int *rev_resultSaIndexLeft, unsigned int *rev_resultSaIndexRight);

void BWTRetrievePositionFromSAIndex(Idx2BWT * idx2BWT, unsigned int saIndex, unsigned int * sequenceId, unsigned int * offset);

int bwt_match_exact(Idx2BWT *bi_bwt, ubyte_t *seq, int len, bwtint_t *sa_begin,
        bwtint_t *sa_end, bwtint_t *rev_sa_begin, bwtint_t *rev_sa_end);

int bwt_extend_exact(const Idx2BWT *bi_bwt, ubyte_t *seq, int len, int *leav_len, int type,
        bwtint_t *sa_begin, bwtint_t *sa_end, bwtint_t *rev_sa_begin, bwtint_t *rev_sa_end);

#endif
