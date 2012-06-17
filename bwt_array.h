#ifndef __BWT_ARRAY_T_
#define __BWT_ARRAY_T_

#include "2BWT-Interface.h"

/*
 * two conjoint bwtint_ts make up an element,
 * first bwtint_t num represent intron position,
 * second bwtint_t num save the hit count and intron pos type.
 * The array is ordered
 * @intron pos type : whether the pos is start of the intron
 * @hit count: how many mapped seq is spliced by the intron
 *
 * an element:
 * ------------------------------------------
 *  intron pos
 * ------------------------------------------
 *  intron pos type| intron type | hit count
 *  1 bit          | 2 bit       | 28 bit
 * ------------------------------------------
 */

typedef struct _array{
    int size; // allocated elements
    int cur_cnt; // used elements
    int last_mid; // last time searched mid num
    bwtint_t *data; // address of elements
} bwt_array_t;

/*
 * init an array, the initial array size is 10000;
 */
inline bwt_array_t *bwt_array_init();

// destory the array
inline void bwt_array_destory(bwt_array_t *arr);

/*
 * the array size is not enough, realloc the array
 */
inline void bwt_array_resize(bwt_array_t *arr);

/*
 * insert _pos into a suitable position, the array is ordered
 */
inline void bwt_array_insert(bwt_array_t *arr, bwtint_t _pos, int intron_type);

/*
 * find a suitable intron position by _pos, direction and seq len
 */
int bwt_find_split_pos_by_record(bwt_array_t *arr, bwtint_t _pos, int len, int is_backward);

#endif
