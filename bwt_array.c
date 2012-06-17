#include "bwt_array.h"

inline bwt_array_t *bwt_array_init()
{
    bwt_array_t *arr;
    arr = (bwt_array_t *)malloc(sizeof(bwt_array_t));
    arr->size = 10000;
    arr->cur_cnt = 0;
    arr->last_mid = 0;
    arr->data = (bwtint_t *) calloc(arr->size * 2 , sizeof(bwtint_t));
    memset(arr->data, 0, arr->size * 2 * sizeof(bwtint_t));
    return arr;
}

inline void bwt_array_destory(bwt_array_t *arr)
{
    free(arr->data);
    free(arr);
}

inline void bwt_array_resize(bwt_array_t *arr)
{
    arr->size <<= 1;
    arr->data = (bwtint_t *)realloc(arr->data, arr->size *2 *sizeof(bwtint_t));
    if(arr->data == NULL) {
        fprintf(stderr, "arr->data can not realloc\n");
        exit(-1);
    }
}

void bwt_array_insert(bwt_array_t *arr, bwtint_t _pos, int intron_type)
{
    int i;
    if(arr->data[arr->last_mid * 2] == _pos)
        arr->data[arr->last_mid * 2 + 1] += 1;
    else{
        if(arr->cur_cnt == arr->size)
            bwt_array_resize(arr);
        if(arr->data[arr->last_mid * 2] > _pos){
            // relocate mid
            while(arr->data[(arr->last_mid - 1) * 2] > _pos)
                arr->last_mid -= 1;
            // last_mid is a insert position
            for(i = arr->cur_cnt; i >= arr->last_mid; --i) {
                arr->data[ (i+1) * 2] = arr->data[i * 2];
                arr->data[(i+1)*2 + 1] = arr->data[i * 2 + 1];
            }
            i += 1;
            arr->data[i * 2] = _pos;
            // _pos sit infront of intron start
            arr->data[i * 2 + 1] = intron_type << 29;
        } else {
            // relocate mid 
            while(arr->data[arr->last_mid + 1 < _pos])
                arr->last_mid += 1;

            for(i = arr->cur_cnt; i > arr->last_mid; --i) {
                arr->data[ (i+1) * 2] = arr->data[i * 2];
                arr->data[(i+1)*2 + 1] = arr->data[i * 2 + 1];
            }
            i += 1;
            arr->data[i * 2] = _pos;
            // _pos sit behand intron end
            arr->data[i * 2 + 1] = 0x80000000 | (intron_type << 29);
        }
    }
}

// result: possible splicing site of seq, relative to _pos
// is_backward : represent that _pos should be sitted in the back of intron
int bwt_find_split_pos_by_record(bwt_array_t *arr, bwtint_t _pos, int len, int is_backward)
{
    // first use binary search, then locate the possible intron
    // it may have more than one possible introns, but only return one
    int low, high, mid;
    if(arr->cur_cnt == 0)
        return -1;

    low = 0;
    high = arr->cur_cnt - 1;
    mid = (low + high)/2;
    if(is_backward == 0){
        if(arr->data[arr->cur_cnt * 2] < _pos)
            return -1;
        while(high > low){
            if(arr->data[mid * 2] > _pos)
                high = mid;
            else if(arr->data[mid * 2] < _pos)
                low = mid + 1;
            mid = (low + high)/2;
        }
        // mid is the first element, whose value is larger than _pos
        arr->last_mid = mid;
        while(_pos >= arr->data[mid * 2] && _pos + len > arr->data[mid * 2])
        {
            // filter the result
            if(arr->data[mid * 2 +1] >> 31 == 0){
                arr->last_mid = mid;
                return arr->data[mid * 2] - _pos;
            }
            ++mid;
        }
        return -1;
    } else {
        if(arr->data[0] > _pos)
            return -1;
        while(high > low){
            if(arr->data[mid * 2] < _pos)
                low = mid;
            else
                high = mid - 1;
            mid = (high + low)/2;
        }
        arr->last_mid = mid;
        // mid is the first element, whose value is less than _pos
        while(_pos >= arr->data[mid * 2] && _pos + len > arr->data[mid * 2])
        {
            if(arr->data[mid * 2 +1] >> 31 == 0){
                arr->last_mid = mid;
                return _pos - arr->data[mid * 2];
            }
            ++mid;
        }
        return -1;
    }
}
