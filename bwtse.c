#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "2BWT-Interface.h"
#include "stdaln.h"
#include "bwtse.h"
#include "bwtaln.h"
//#include "bntseq.h"
#include "utils.h"
#include "kstring.h"

int g_log_n[256];
char *bwt_rg_line, *bwt_rg_id;

int max_tot_cnt = 30;
bwt_pos_t *pos_array;

void bwt_print_sam_PG();
void bwt_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi)
{
    int i, cnt, best;
    if(n_aln == 0){
        s->type = BWA_TYPE_NO_MATCH;
        s->c1 = s->c2 = 0;
        return;
    }

    //TODO when the aln is splicing mapping?

    if(set_main){
        best = aln[0].score;
        for(i = cnt = 0; i < n_aln; ++i){
            const bwt_aln1_t *p = aln + i;
            if(p->score > best)
                break;
            if(drand48()*(p->l - p->k + 1 + cnt) > (double) cnt){
                s->n_mm = p->n_mm;
                s->n_gapo = p->n_gapo;
                s->n_gape = p->n_gape;
                s->score = p->score;
                s->start = p->start;
                s->end = p->end;
                s->sa = p->k + (bwtint_t)((p->l - p->k + 1) * drand48());
                s->strand = p->strand;
            }
            cnt += p->l - p->k + 1;
        }
        s->c1 = cnt;
        for(; i < n_aln; ++i)
            cnt += aln[i].l - aln[i].k + 1;
        s->c2 = cnt - s->c1;
        s->type = s->c1 > 1? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
    }
    if(n_multi){
        int k, rest, n_occ, z = 0;
        for(k = n_occ = 0; k < n_aln; ++k){
            const bwt_aln1_t *q = aln + k;
            n_occ += q->l - q->k + 1;
        }
        if(s->multi)
            free(s->multi);
        if(n_occ > n_multi + 1){
            //there are too many hits, return none
            s->multi = 0;
            s->n_multi = 0;
            return ;
        }

        //find one additional for -> sa
        rest = n_occ > n_multi + 1? n_multi + 1 : n_occ;
        s->multi = calloc(rest, sizeof(bwt_multi1_t));
        for(k = 0; k < n_aln; ++k){
            const bwt_aln1_t *q = aln + k;
            if(q->l - q->k + 1 <= rest){
                bwtint_t l;
                for(l = q->k; l <= q->l; ++l){
                    s->multi[z].start = q->start;
                    s->multi[z].end = q->end;
                    s->multi[z].strand = q->strand;
                    s->multi[z].pos = l;
                    s->multi[z].gap = q->n_gapo + q->n_gape;
                    s->multi[z++].mm = q->n_mm;
                }
                rest -= q->l - q->k + 1;
            } else {
                //Random sample
                int j, i;
                for(j = rest, i = q->l - q->k + 1, k = 0; j > 0; --j){
                    double p = 1.0, x = drand48();
                    while(x < p)
                        p -= p*j/(i--);
                    s->multi[z].start = q->start;
                    s->multi[z].end = q->end;
                    s->multi[z].strand = q->strand;
                    s->multi[z].pos = q->l - i;
                    s->multi[z].gap = q->n_gapo + q->n_gape;
                    s->multi[z++].mm = q->n_mm;
                }
                rest = 0;
                break;
            }
        }
        s->n_multi = z;
    }
}
// TODO if alignment is splicing mapping, retrieve every aln,
// get seq id and offset, select the most suitable combination 

void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s)
{
    bwt_aln2seq_core(n_aln, aln, s, 1, 0);
}

int bwa_approx_mapQ(const bwa_seq_t *p, int mm)
{
    int n;
    if(p->c1 == 0) return 23;
    if(p->c1 > 1) return 0;
    if(p->n_mm == mm) return 25;
    if(p->c2 == 0) return 37;
    n = (p->c2 >= 255)? 255 : p->c2;
    return (23 < g_log_n[n])? 0 : 23 - g_log_n[n];
}

/*
 * Derive the actual position in the read from the given suffix array
 * coordinates. Note that the position will be approximate based on 
 * whether indels appear in the read and whether calculations are
 * performed from the start or end of the read.
 */
void bwa_cal_pac_pos_core(const Idx2BWT *bi_bwt, bwa_seq_t *seq, const int max_mm, const float fnr)
{
    int max_diff;
    if(seq->type != BWA_TYPE_UNIQUE && seq->type != BWA_TYPE_REPEAT)
        return ;
    max_diff = fnr > 0.0 ? bwa_cal_maxdiff(seq->len, BWA_AVG_ERR, fnr) : max_mm;
    seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff);
    BWTRetrievePositionFromSAIndex(bi_bwt, seq->sa, &seq->seq_id, &seq->pos);
    seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff);
}

// use greedy algorithm, pos_res is the result
// 可能不是最优解
int bwt_find_approx_aln(bwt_pos_t *pos_array, bwt_pos_t *pos_res, int n_aln)
{
    int i, j,k;
    int index = 0;
    bwt_pos_t *tmp_arr = (bwt_pos_t *)calloc(n_aln, sizeof(bwt_pos_t));
    int flag, stat = 0;
    j = 0;
    bwtint_t len_intron;
    bwtint_t min_len_intron;
    while(j < n_aln && (pos_array+j)->r_aln == 0){
        index = j;
        memcpy(tmp_arr, pos_array+j, sizeof(bwt_pos_t));
        flag = 0;
        for(k = 1; k < n_aln; ++k){
            while((pos_array+index)->r_aln != k){
                ++index;
                continue;
            }
            int k_aln_start = index;
            while((pos_array+index)->r_aln == k){
                // TODO 正负链的坐标是不是有所不同
                if(k_aln_start == index || (((pos_array+index)->seq_id == (tmp_arr+k-1)->seq_id) &&
                        ((pos_array+index)->strand == (tmp_arr+k-1)->strand) &&
                        ((pos_array+index)->pos > (tmp_arr+k-1)->pos) && 
                        ((pos_array+index)->pos < (tmp_arr+k)->pos))){
                    memcpy(tmp_arr+k, pos_array + index, sizeof(bwt_pos_t));
                    ++index;
                    flag = 1;
                } else
                    ++index;
            }
            // 在第k个aln 没有符合标准的
            if(0 == flag)
                break;
        }
        //没有符合条件的
        if(flag == 0)
            continue;
        //cal total length of intron
        len_intron = 0;
        for(i = 1; i != n_aln; ++i)
            len_intron += (tmp_arr+i)->pos - ((tmp_arr+i-1)->pos + 
                    (tmp_arr+i-1)->end-(tmp_arr+i-1)->start + 1);
        if(stat == 0 || len_intron < min_len_intron){
            memcpy(pos_res, tmp_arr, sizeof(bwt_pos_t)*n_aln);
            min_len_intron = len_intron;
            stat = 1;
        }
        ++j;
    }
    free(tmp_arr);
    return stat;
}

// by bwt_multi1_t struct, multi save aprart of pos info of splicing mapping
void bwt_aln2pos_splicing(Idx2BWT *bi_bwt, bwa_seq_t *seq, int n_aln, bwt_aln1_t *aln,
        int max_diff, float fnr)
{
    int i, k, cnt = 0;
    int index = 0;
    int strand;
    int nm, res_flag;
    unsigned int j, seq_id;
    bwt_pos_t *pos_info;
    bwt_aln1_t *tmp_aln;
    tmp_aln = aln;
    memset(pos_array, 0, sizeof(bwt_pos_t) * max_tot_cnt);
    for(i = 0; i != n_aln; ++i){
        //TODO if cnt to large, discard id
        cnt += tmp_aln->l - tmp_aln->k + 1;
        nm = tmp_aln->n_mm + tmp_aln->n_gapo + tmp_aln->n_gape;
        if(cnt > max_tot_cnt){
            max_tot_cnt = cnt * 2;
            pos_array = (bwt_pos_t *)realloc(pos_array, sizeof(bwt_pos_t) * max_tot_cnt);
        }
        for(j = tmp_aln->k; j <= tmp_aln->l; ++j){
            BWTRetrievePositionFromSAIndex(bi_bwt, j, &seq_id, &(pos_array+index)->pos);
            (pos_array+index)->strand = tmp_aln->strand;
            (pos_array+index)->seq_id = seq_id;
            pos_array[index].start = tmp_aln->start;
            pos_array[index].end = tmp_aln->end;
            pos_array[index].r_aln = i;
            pos_array[index].nm = nm;
            index += 1;
        }
        ++tmp_aln;
    }
    //find the closest pos to form the alignment
    int r_aln;
    pos_info = (bwt_pos_t *)calloc(n_aln, sizeof(bwt_pos_t));
    if(n_aln == 1) {
        memcpy(pos_info, pos_array, sizeof(bwt_pos_t));
        res_flag = 1;
    } else
        res_flag = bwt_find_approx_aln(pos_array, pos_info, n_aln);
    if(res_flag == 1){
        seq->strand = pos_info->strand;
        seq->seq_id = pos_info->seq_id;
        seq->pos = pos_info->pos;
        seq->start = pos_info->start;
        seq->end = pos_info->end;
        seq->n_mm = aln->n_mm;
        seq->n_gapo = aln->n_gapo;
        seq->n_gape = aln->n_gape;
        seq->n_aln = 1;
        bwt_multi1_t *multi = (bwt_multi1_t *)calloc(n_aln-1, sizeof(bwt_multi1_t));
        seq->multi = multi;
        seq->n_multi = n_aln -1;
        for(i = 1; i < n_aln; ++i){
            r_aln = (pos_info + i)->r_aln;
            multi->strand = (pos_info + i)->strand;
            multi->seq_id = (pos_info + i)->seq_id;
            multi->pos = (pos_info + i)->pos;
            multi->start = (pos_info + i)->start;
            multi->end = (pos_info + i)->end;
            multi->mm = (aln+r_aln)->n_mm;
            multi->gap = (aln+r_aln)->n_gapo + (aln+r_aln)->n_gape;
        }
    } else {
        seq->type = BWA_TYPE_NO_MATCH;
        seq->c1 = seq->c2 = 0;
    }
    free(pos_info);
}
void bwa_cal_pac_pos(const Idx2BWT *bi_bwt, int n_seqs, bwa_seq_t *seq, int max_mm, float fnr)
{
    int i, j, n_multi;
    for(i = 0; i != n_seqs; ++i){
        bwa_seq_t *p = &seq[i];
        // do not deal with splicing mapping
        if(p->type == BWA_TYPE_SPLICING)
            continue;
        bwa_cal_pac_pos_core(bi_bwt, p, max_mm, fnr); 
        for(j = n_multi = 0; j < p->n_multi; ++j){
            bwt_multi1_t *q = p->multi + j;
            BWTRetrievePositionFromSAIndex(bi_bwt, q->pos, &q->seq_id, &q->pos);
            if(q->pos != p->pos)
                p->multi[n_multi++] = *q;
        }
        p->n_multi = n_multi;
    }
}

/*
 * is_end_correct == 1 if (*pos + len)gives the correct coordinate on
 * forward strand. This happens when p->pos is calculated by
 * bwa_cal_pac_pos(). is_end_correct == 0 if(*pos) gives the correct 
 * coordinate. This happens only for color-converted alignment.
 */

// TODO how about splicing mapping, cal cigar use N flag
// len will return leaving length, which without deletion on two end
static bwa_cigar_t *refine_gapped_core(const HSP *hsp, int len, const ubyte_t *seq, bwtint_t *_pos,
        int ext, int *n_cigar, int is_end_correct, int *start, int *end)
{
    bwa_cigar_t *cigar = 0;
    ubyte_t *ref_seq;
    int l = 0, path_len, ref_len;
    AlnParam ap = aln_param_bwa;
    path_t *path;
    int64_t k, __pos = *_pos;

    //ext : the num of gaps
    ref_len = len + abs(ext);
    if(ext > 0){
        ref_seq = (ubyte_t *)calloc(ref_len, 1);
        for(k = __pos; k < __pos + ref_len && k < hsp->dnaLength; ++k)
            //TODO what does packedDNA saved?
            ref_seq[l++] = hsp->packedDNA[k>>4] >> ((~k&15)<<1) & 3;
    } else{
        int64_t x = __pos + (is_end_correct? len : ref_len);
        ref_seq = (ubyte_t *)calloc(ref_len, 1);
        for(l = 0, k = x - ref_len > 0? x - ref_len : 0; k < x && k < hsp->dnaLength; ++k)
            ref_seq[l++] = hsp->packedDNA[k>>4] >> ((~k&15)<<1) &3;
    }
    path = (path_t *)calloc(l + (len), sizeof(path_t));

    aln_global_core(ref_seq, l, (ubyte_t *)seq, len, &ap, path, &path_len);
    cigar = bwa_aln_path2cigar(path, path_len, n_cigar);

    if(ext < 0 && is_end_correct){
        //fix coordinate for reads mapped to the forward strand
        for(l = k = 0; k < *n_cigar; ++k){
            if(__cigar_op(cigar[k]) == FROM_D)
                l -= __cigar_len(cigar[k]);
            else if(__cigar_op(cigar[k]) == FROM_I)
                l += __cigar_len(cigar[k]);
        }
        __pos += l;
    }

    if(__cigar_op(cigar[0]) == FROM_D){
        //deletion at 5'- end
        __pos += __cigar_len(cigar[0]);
        *start += __cigar_len(cigar[0]);
        for(k = 0; k < *n_cigar - 1; ++k)
            cigar[k] = cigar[k+1];
        --(*n_cigar);
    }
    if(__cigar_op(cigar[*n_cigar - 1]) == FROM_D){
        //deletion at the 3'end
        *end -= __cigar_len(cigar[*n_cigar - 1]);
        --(*n_cigar);
    }

    //change "I" at either end of the read to S. just in case. This should rarely happen.
    // TODO if it happened, does it have impact on start and end;
    if(__cigar_op(cigar[*n_cigar - 1]) == FROM_I)
            cigar[*n_cigar - 1] = __cigar_create(FROM_S, (__cigar_len(cigar[*n_cigar-1])));

    if(__cigar_op(cigar[0]) == FROM_I)
            cigar[0] = __cigar_create(FROM_S, (__cigar_len(cigar[0])));

    *_pos = (bwtint_t)__pos;
    free(ref_seq);
    free(path);
    return cigar;
}

char *bwa_cal_md1(HSP *hsp, int n_cigar, bwa_cigar_t *cigar, int len, bwtint_t pos,
         ubyte_t *seq, kstring_t *str, int *_nm)
{
    bwtint_t x, y;
    int z, u, c, nm = 0;
    str->l = 0; //reset 
    x = pos; y = 0;
    if(cigar){
        int k, l;
        for(k = u = 0; k < n_cigar; ++k){
            l = __cigar_len(cigar[k]);
            if(__cigar_op(cigar[k]) == FROM_M){
                for(z = 0; z < l && x+z < hsp->dnaLength; ++z){
                    c = hsp->packedDNA[(x+z)>>4] >> ((~(x+z)&15)<<1) & 3;
                    if(c > 3|| seq[y+z] > 3 || c != seq[y+z]){
                        ksprintf(str, "%d", u);
                        kputc("ACGTN"[c], str);
                        ++nm;
                        u = 0;
                    } else
                        ++u;
                } 
                x += l; y += l;
            } else if(__cigar_op(cigar[k]) == FROM_I || __cigar_op(cigar[k]) == FROM_S){
                y += l;
                if(__cigar_op(cigar[k]) == FROM_I)
                    nm += l;
            }else if(__cigar_op(cigar[k]) == FROM_D){
                ksprintf(str, "%d", u);
                kputc('^', str);
                for(z = 0; z < l && x+z < hsp->dnaLength; ++z)
                    kputc("ACGT"[hsp->packedDNA[(x+z)>>4]>>((~(x+z)&15)<<1) & 3], str);
                u = 0;
                x += l; nm += l;
            }
        }
    } else {
        //no gaps
        for(z = u = 0; z < (bwtint_t)len; ++z){
            c = hsp->packedDNA[(x+z)>>4] >> ((~(x+z)&15)<<1) & 3;
            if(c >3 || seq[y+z] > 3|| c != seq[y+z]){
                ksprintf(str, "%d", u);
                kputc("ACGTN"[c], str);
                ++nm;
                u = 0;
            }else 
                ++u;
        }
    }
    ksprintf(str, "%d", u);
    *_nm = nm;
    return strdup(str->s);
}

void bwa_correct_trimmed(bwa_seq_t *s)
{
    if(s->len == s->full_len)
        return;
    if(s->strand == 0){
        //forward
        if(s->cigar && __cigar_op(s->cigar[s->n_cigar-1]) == FROM_S){
            //the last is S
            s->cigar[s->n_cigar-1] += s->full_len - s->len;
        } else {
            if(s->cigar == 0){
                s->n_cigar = 2;
                s->cigar = calloc(s->n_cigar, sizeof(bwa_cigar_t));
                s->cigar[0] = __cigar_create(FROM_M, s->len);
            } else {
                ++s->n_cigar;
                s->cigar = realloc(s->cigar, s->n_cigar * sizeof(bwa_cigar_t));
            }
            s->cigar[s->n_cigar - 1] = __cigar_create(FROM_S, (s->full_len - s->len));
        }
    } else{
        //reverse strand
        if(s->cigar && __cigar_op(s->cigar[0]) == FROM_S)
            s->cigar[0] += s->full_len - s->len;
        else{
            if(s->cigar == 0){
                s->n_cigar = 2;
                s->cigar = calloc(s->n_cigar, sizeof(bwa_cigar_t));
                s->cigar[1] = __cigar_create(FROM_M, s->len);
            }else{
                ++s->n_cigar;
                s->cigar = realloc(s->cigar, s->n_cigar * sizeof(bwa_cigar_t));
                memmove(s->cigar + 1, s->cigar, (s->n_cigar - 1) *sizeof(bwa_cigar_t));
            }
            s->cigar[0] = __cigar_create(FROM_S, (s->full_len - s->len));
        }
    }
    s->len = s->full_len;
}

void bwa_refine_gapped(const HSP *hsp, int n_seqs, bwa_seq_t *seqs, HSP *hsp_cs)
{
    int i, j, n_cigar;
    kstring_t *str;
    ubyte_t *str_tmp;
    int max_len = 100;
    str_tmp = (ubyte_t *)calloc(max_len, sizeof(ubyte_t));

    for(i = 0; i != n_seqs; ++i){
        bwa_seq_t *s = seqs + i;
        //TODO Does it need reverse
        //seq_reverse(s->len, s->seq, 0);
        if(s->type != BWA_TYPE_SPLICING){
            // cal cigar for common mapping
            for(j = 0; j < s->n_multi; ++j){
                bwt_multi1_t *q = s->multi + j;
                int n_cigar;
                // because of it does not need cal cigar 
                if(q->gap == 0)
                    continue;
                q->cigar = refine_gapped_core(hsp, s->len, q->strand? s->rseq:s->seq, &q->pos,
                        (q->strand? 1 : -1)*q->gap, &n_cigar, 1, &q->start, &q->end);
                q->n_cigar = n_cigar;
            }
            if(s->type == BWA_TYPE_NO_MATCH || s->type == BWA_TYPE_MATESW || s->n_gapo == 0)
                continue;
            s->cigar = refine_gapped_core(hsp, s->len, s->strand? s->rseq: s->seq, &s->pos,
                    (s->strand? 1: -1)*(s->n_gapo + s->n_gape), &s->n_cigar, 1, &s->start, &s->end);
        } else{
            //for splicing mapping
            int tot_cigar = 0;
            strncpy(str_tmp, (s->strand? s->rseq:s->seq) + s->start, s->end - s->start + 1 );
            s->cigar = refine_gapped_core(hsp, s->end - s->start+1, str_tmp, &s->pos,
                    (s->strand? 1:-1)*(s->n_gapo + s->n_gape), &s->n_cigar, 1, &s->start, &s->end);
            tot_cigar = s->n_cigar;
            for(j = 0; j < s->n_multi; ++j){
                bwt_multi1_t *q = s->multi + j;
                int n_cigar;
                strncpy(str_tmp, (s->strand? s->rseq:s->seq) + q->start, q->end-q->start+1);
                q->cigar = refine_gapped_core(hsp, q->end-q->start+1, str_tmp, &q->pos,
                        (s->strand? 1:-1)*q->gap, &n_cigar, 1, &q->start, &q->end);
                q->n_cigar = n_cigar;
                tot_cigar += q->n_cigar;
            }
            // merge aln to splicing type
            tot_cigar += 2;
            int start, end;
            s->cigar = (bwa_cigar_t *)realloc(s->cigar, sizeof(bwa_cigar_t)*(tot_cigar + s->n_multi));
            start = s->start;
            end = s->end;
            // m_cigar = s->n_cigar;
            for(j = 0; j < s->n_multi; ++j){
                bwt_multi1_t *tmp = s->multi + j;
                int ov_len = end - tmp->start;
                // remove overlapping part in alignment
                if(ov_len > 0){
                    if(__cigar_len(tmp->cigar[0]) > ov_len){
                        __cigar_mod_len(tmp->cigar[0], __cigar_len(tmp->cigar[0]) -ov_len);
                        tmp->start += ov_len;
                        tmp->pos += ov_len;
                    }
                    else if(__cigar_len(s->cigar[s->n_cigar-1]) > ov_len){
                        __cigar_mod_len(s->cigar[s->n_cigar-1], 
                                __cigar_len(s->cigar[s->n_cigar-1] - ov_len));
                        s->end -= ov_len;
                    } else if(__cigar_len(tmp->cigar[0] > ov_len>>1 &&
                                __cigar_len(s->cigar[s->n_cigar-1]) > ov_len>>1)){
                        __cigar_mod_len(s->cigar[s->n_cigar-1], __cigar_len(s->cigar[s->n_cigar-1])-(ov_len>>1));
                        __cigar_mod_len(tmp->cigar[0], __cigar_len(tmp->cigar[0])-ov_len+(ov_len>>1));
                        s->end -= (ov_len>>1);
                        tmp->start += ov_len - (ov_len>>1);
                        tmp->pos += ov_len - (ov_len>>1);
                    }else
                        // too much overlapping part, discard it
                        continue;
                    s->cigar[s->n_cigar]= __cigar_create(FROM_N, tmp->pos - s->pos - s->end);
                    memcpy(&s->cigar[s->n_cigar+1], tmp->cigar, sizeof(bwa_cigar_t)*tmp->n_cigar);
                    s->n_cigar += tmp->n_cigar + 1;
                } else if(ov_len < 0){
                    // having unmapped part in read, handle it like insert
                    s->cigar[s->n_cigar]= __cigar_create(FROM_I, tmp->start - s->end);
                    s->cigar[s->n_cigar + 1]= __cigar_create(FROM_I, tmp->pos - s->pos - tmp->start);
                    memcpy(&s->cigar[s->n_cigar+2], tmp->cigar, sizeof(bwa_cigar_t)*tmp->n_cigar);
                    s->n_cigar += tmp->n_cigar + 2;
                }
            }
            s->n_multi=0;
            free(s->multi);
        }
    }
/*     if(hsp_cs) {
 *         // TODO splicing mapping is not support color space currently
 *         // in color space
 *         for(i = 0; i < n_seqs; ++i){
 *             bwa_seq_t *s = seqs + i;
 *             //TODO
 *             bwa_cs2nt_core(s, hsp->dnaLength, hsp_cs);
 *             for(j = 0; j < s->n_multi; ++j){
 *                 bwt_multi1_t *q = s->multi + j;
 *                 int n_cigar;
 *                 if(q->gap == 0)
 *                     continue;
 *                 free(q->cigar);
 *                 q->cigar = refine_gapped_core(hsp_cs, s->len, q->strand? s->rseq: s->seq, &q->pos,
 *                         (q->strand? 1:-1) * q->gap, &n_cigar, 0,&q->start, &q->end);
 *                 q->n_cigar = n_cigar;
 *             }
 *             if(s->type != BWA_TYPE_NO_MATCH || s->cigar){
 *                 //update cigar
 *                 free(s->cigar);
 *                 s->cigar = refine_gapped_core(hsp_cs, s->len, s->strand? s->rseq: s->seq, &s->pos,
 *                         (s->strand?1:-1) * (s->n_gapo + s->n_gape), &s->n_cigar, 0, &s->start, &s->end);
 *             }
 *         }
 *     }
 */

    //generate MD tag
    str = (kstring_t *)calloc(1, sizeof(kstring_t));
    for(i = 0; i != n_seqs; ++i){
        bwa_seq_t *s = seqs + i;
        if(s->type != BWA_TYPE_NO_MATCH && s->type != BWA_TYPE_SPLICING){
            int nm;
/*             s->md = bwa_cal_md1(hsp_cs? hsp_cs: hsp, s->n_cigar, s->cigar, s->len, s->pos, s->strand? s->rseq: s->seq,
 *                      str, &nm);
 */
            s->md = bwa_cal_md1(hsp, s->n_cigar, s->cigar, s->len, s->pos, s->strand? s->rseq: s->seq,
                     str, &nm);
            s->nm = nm;
        }
    }
    free(str->s);
    free(str);

    //correct for trimmed reads
/*     if(!hsp_cs){
 *         //It is only for Illumina reads
 *         for(i = 0; i < n_seqs; ++i)
 *             bwa_correct_trimmed(seqs + i);
 *     }
 */
    for(i = 0; i < n_seqs; ++i)
        bwa_correct_trimmed(seqs + i);
}

int64_t pos_end(const bwa_seq_t *p)
{
    if(p->cigar){
        int j;
        int64_t x = p->pos;
        for(j = 0; j != p->n_cigar; ++j){
            int op = __cigar_op(p->cigar[j]);
            if(op == 0 || op == 2)
                x += __cigar_op(p->cigar[j]);
        }
        return x;
    } else 
        return p->pos + p->len;
}

int64_t pos_end_multi(const bwt_multi1_t *p, int len)
{
    if(p->cigar){
        int j;
        int64_t x = p->pos;
        for(j = 0; j != p->n_cigar; ++j){
            int op = __cigar_op(p->cigar[j]);
            if(op == 0 || op == 2)
                x += __cigar_len(p->cigar[j]);
        }
        return x;
    } else
        return p->pos + len;
}
static int64_t pos_5(const bwa_seq_t *p)
{
    if(p->type != BWA_TYPE_NO_MATCH)
        return p->strand? pos_end(p) : p->pos;
    return -1;
}

void bwa_print_sam1(const HSP *hsp, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2)
{
    int seq_id, j;
    if(p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)){
        int nn, am = 0, flag = p->extra_flag;
        char XT;

        if(p->type == BWA_TYPE_NO_MATCH){
            p->pos = mate->pos;
            p->strand  = mate->strand;
            flag |= SAM_FSU;
            j = 1;
        } else 
            j = pos_end(p) - p->pos;    //j is the length of the reference in the alignment

        //TODO cal nn
        seq_id = p->seq_id;
        if(p->type != BWA_TYPE_NO_MATCH && p->pos + j - hsp->ambiguity[seq_id].startPos > hsp->ambiguity[seq_id].rightOfEndPos)
            flag |= SAM_FSU; //flag UNMAP as this alignment bridges two adjacent reference sequences
        //update flag and print it
        if(p->strand)
            flag |= SAM_FSR;
        if(mate){
            if(mate->type != BWA_TYPE_NO_MATCH){
                if(mate->strand)
                    flag |= SAM_FMR;
            } else 
                flag |= SAM_FMU;
        }
        err_printf("%s\t%d\t%s\t", p->name, flag, hsp->annotation[seq_id].text);
        err_printf("%d\t%d\t",(int)(p->pos - hsp->seqOffset[seq_id].startPos + 1), p->mapQ);

        //print CIGAR
        if(p->cigar){
            for(j = 0; j != p->n_cigar; ++j)
                err_printf("%d%c", __cigar_len(p->cigar[j]), "MIDNSHP=X"[__cigar_op(p->cigar[j])]);
        } else if(p->type == BWA_TYPE_NO_MATCH)
            err_printf("*");
        else
            err_printf("%dM", p->len);

        //print mate coordinate
        if(mate && mate->type != BWA_TYPE_NO_MATCH){
            int m_seqid, m_is_N;
            long long isize;
            am = mate->seQ < p->seQ? mate->seQ : p->seQ;    //smaller single-end mapping quality
            //redundant calculation here, but should not matter too much
            //TODO cal nn, 
            //BWTRetrievePositionFromSAIndex(bi_bwt, mate->pos, &mate->strand, &mate->seq_id, &mate->pos);
            err_printf("\t%s\t",(seq_id == mate->seq_id)?"=":hsp->annotation[mate->seq_id].text);
            isize = (seq_id == mate->seq_id)? pos_5(mate) - pos_5(p) : 0;
            if(p->type == BWA_TYPE_NO_MATCH)
                isize = 0;
            err_printf("%d\t%lld\t", (int)(mate->pos - hsp->seqOffset[mate->seq_id].startPos + 1), isize);
        }else if(mate)
            err_printf("\t=\t%d\t0\t",(int)(p->pos - hsp->seqOffset[seq_id].startPos + 1));
        else
            err_printf("\t*\t0\t0\t");

        //print sequence and quality
        if(p->strand == 0)
            for(j = 0; j != p->full_len; ++j)
                putchar("ACGTN"[(int)p->seq[j]]);
        else
            for(j = 0; j != p->full_len; ++j)
                putchar("TGCAN"[p->seq[p->full_len - 1 -j]]);
        putchar('\t');
        
        if(p->qual){
            if(p->strand)
                seq_reverse(p->len, p->qual,0); //reverse quality
            err_printf("%s", p->qual);
        }else
            printf("*");

        if(bwt_rg_id)
            err_printf("\tRG:Z:%s", bwt_rg_id);
        if(p->bc[0])
            err_printf("\tBC:Z:%s",p->bc);
        if(p->clip_len < p->full_len)
            err_printf("\tXC:i:%d", p->clip_len);
        
        if(p->type != BWA_TYPE_NO_MATCH){
            int i;
            //calculate XT tag
            XT = "NURM"[p->type];
            if(nn > 10)
                XT = 'N';
            //print tags
            err_printf("\tXT:A:%c\t%s:i:%d", XT,(mode & BWA_MODE_COMPREAD)?"NM":"CM",p->nm);
            if(nn)
                err_printf("\tXN:i:%d", nn);
            if(mate)
                err_printf("\tSM:i:%d\tAM:i:%d", p->seQ, am);
            if(p->type != BWA_TYPE_MATESW){
                //X0 and X1 are not available for this type of alignment
                err_printf("\tX0:i:%d", p->c1);
                if(p->c1 <= max_top2)
                    err_printf("\tX1:i:%d", p->c2);
            }
            err_printf("\tXM:i:%d\tXO:i:%d\tXG:i:%d", p->n_mm, p->n_gapo, p->n_gapo+p->n_gape);
            if(p->md)
                err_printf("\tMD:Z:%s", p->md);
            //print multiple hits
            if(p->n_multi){
                err_printf("\tXA:Z:");
                for(i = 0; i < p->n_multi; ++i){
                    bwt_multi1_t *q = p->multi + i;
                    int k;
                    j = pos_end_multi(q, p->len) - q->pos;
                    //TODO when multiple hits, get pos before print to sam files
                    err_printf("%s,%c%d",hsp->annotation[q->seq_id].text, q->strand? '-' : '+',
                            (int)(q->pos - hsp->seqOffset[q->seq_id].startPos +1));
                    if(q->cigar){
                        for(k = 0; k < q->n_cigar; ++k)
                            err_printf("%d%c", __cigar_len(q->cigar[k]), "MIDS"[__cigar_op(q->cigar[k])]);
                    } else 
                        err_printf(",%d;", q->gap + q->mm);
                }
            }
        }
        putchar('\n');
    } else {
        //this read has no match
        ubyte_t *s = p->strand? p->rseq : p->seq;
        int flag = p->extra_flag | SAM_FSU;
        if(mate && mate->type == BWA_TYPE_NO_MATCH)
            flag |= SAM_FMU;
        err_printf("%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", p->name, flag);
        for(j = 0; j != p->len; ++j)
            putchar("ACGTN"[(int)s[j]]);
        putchar('\t');
        if(p->qual){
            if(p->strand)
                seq_reverse(p->len, p->qual, 0);    //reverse quality
            err_printf("%s", p->qual);
        } else
            err_printf("*");
        if(bwt_rg_id)
            err_printf("\tRG:Z:%s", bwt_rg_id);
        if(p->bc[0])
            err_printf("\tBC:Z:%s", p->bc);
        if(p->clip_len < p->full_len)
            err_printf("\tXC:i:%d", p->clip_len);
        putchar('\n');
    }
}

void bwt_print_sam_SQ(const HSP *hsp)
{
    int i;
    for(i = 0; i < hsp->numOfSeq; ++i)
        err_printf("@SQ\tSN:%s\tLN:%d\n", hsp->annotation[i].text,
                hsp->seqOffset[i].endPos - hsp->seqOffset[i].startPos + 1);
    if(bwt_rg_line)
        err_printf("%s\n", bwt_rg_line);
}

void bwase_initialize()
{
    int i;
    for(i = 1; i != 256; ++i)
        g_log_n[i] = (int)(4.343 * log(i) + 0.5);
}

char *bwa_escape(char *s)
{
    char *p, *q;
    for( p = q = s; *p; ++p){
        if(*p == '\\'){
            ++p;
            if(*p == 't')       *q++ = '\t';
            else if(*p == 'n')  *q++ = '\n';
            else if(*p == 'r')  *q++ = '\r';
            else if(*p == '\\') *q++ = '\\';
        } else 
            *q++ = *p;
    }
    *q = '\0';
    return s;
}

int bwa_set_rg(const char *s)
{
    char *p, *q, *r;
    if(strstr(s, "@RG") != s)
        return -1;
    if(bwt_rg_line)
        free(bwt_rg_line);
    if(bwt_rg_id)
        free(bwt_rg_id);
    bwt_rg_line = strdup(s);
    bwt_rg_id = 0;
    bwa_escape(bwt_rg_line);
    p = strstr(bwt_rg_line, "\tID:");
    if(p == 0)
        return -1;
    p += 4;
    for(q = p; *q && *q != '\t' && *q != '\n'; ++q);
    bwt_rg_id = calloc(q - p + 1, 1);
    for(q = p, r = bwt_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
        *r++ = *q;
    return 0;
}

void bwt_sai2sam_se_core(const char *prefix, const char *fn_sa, const char *fn_fa, int n_occ)
{
    extern bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
    int i, n_seqs, tot_seqs = 0, m_aln;
    bwt_aln1_t *aln = 0;
    bwa_seq_t *seqs;
    bwa_seqio_t *ks;
    Idx2BWT *bi_bwt, *bi_bwt_cs;
    clock_t t;
    FILE *fp_sa;
    gap_opt_t opt;

    //initialization
    bwase_initialize();
    // load 2BWT
    char *str = (char *)calloc(strlen(prefix) + 10, 1);
    strcpy(str, prefix);
    strcat(str, ".index");
    bi_bwt = BWTLoad2BWT(str, ".sa");
    free(str);

    srand48(bi_bwt->hsp->randomSeed);
    fp_sa = xopen(fn_sa, "r");

    m_aln = 0;
    bi_bwt_cs = NULL;
    fread(&opt, sizeof(gap_opt_t), 1, fp_sa);
    if(!(opt.mode & BWA_MODE_COMPREAD)){
        // in color space; initialize bidirection bwt
        // load 2BWT
        char *str = (char *)calloc(strlen(prefix) + 20, 1);
        strcpy(str, prefix);
        strcat(str, ".nt.index");
        bi_bwt_cs = BWTLoad2BWT(str, ".sa");
        free(str);
    }

    bwt_print_sam_SQ(bi_bwt->hsp);
    bwt_print_sam_PG();
    //set ks
    ks = bwa_open_reads(opt.mode, fn_fa);
    //init pos_array
    pos_array = (bwt_pos_t *)calloc(max_tot_cnt, sizeof(bwt_pos_t));
    //core loop
    while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt.mode, opt.trim_qual)) != 0){
        tot_seqs += n_seqs;
        t = clock();

        //read alignment
        for(i = 0; i < n_seqs; ++i){
            bwa_seq_t *p = seqs + i;
            int n_aln;
            fread(&n_aln, sizeof(int), 1, fp_sa);
            if(n_aln > m_aln){
                m_aln = n_aln;
                aln = (bwt_aln1_t *)realloc(aln, sizeof(bwt_aln1_t) * m_aln);
            }
            fread(aln, sizeof(bwt_aln1_t), n_aln, fp_sa);
            if(aln->type == BWA_TYPE_SPLICING && n_aln != 0){
                // if the result of splicing mapping
                bwt_aln2pos_splicing(bi_bwt, p, n_aln, aln, opt.max_diff, opt.fnr);
                p->type = BWA_TYPE_SPLICING;
            }
            else
                bwt_aln2seq_core(n_aln, aln, p, 1, n_occ);
        }

        fprintf(stderr, "[bwa_aln_core] convert to sequence coordinate... ");
        bwa_cal_pac_pos(bi_bwt, n_seqs, seqs, opt.max_diff, opt.fnr);
        fprintf(stderr,"%.2f sec\n", (float)(clock() - t)/CLOCKS_PER_SEC);
        t = clock();

        fprintf(stderr, "[bwa_aln_core] refine gapped alignment... ");
        bwa_refine_gapped(bi_bwt->hsp, n_seqs, seqs, bi_bwt_cs == NULL? NULL : bi_bwt_cs->hsp);
        fprintf(stderr,"%.2f sec\n", (float)(clock() - t)/CLOCKS_PER_SEC);
        t = clock();

        fprintf(stderr, "[bwa_aln_core] print alignment... ");
        for(i = 0; i < n_seqs; ++i)
            bwa_print_sam1(bi_bwt->hsp, seqs + i, 0, opt.mode, opt.max_top2);
        fprintf(stderr,"%.2f sec\n", (float)(clock() - t)/CLOCKS_PER_SEC);
        t = clock();

        bwa_free_read_seq(n_seqs, seqs);
        fprintf(stderr,"[bwa_aln_core] %d sequences have been processed.\n", tot_seqs);
    }

    //destory bwt
    bwa_seq_close(ks);
    if(bi_bwt_cs)
        BWTFree2BWT(bi_bwt_cs);
    BWTFree2BWT(bi_bwt);
    fclose(fp_sa);
    free(aln);
}

int bwa_sai2sam_se(int argc, char *argv[])
{
    int c, n_occ = 3;
    while((c = getopt(argc, argv, "hn:f:r:")) >= 0){
        switch(c){
            case 'h': break;
            case 'r':
                if(bwa_set_rg(optarg) < 0){
                    fprintf(stderr, "[%s] malformated @RG line\n", __func__);
                    return 1;
                }
                break;
            case 'n': n_occ = atoi(optarg); break;
            case 'f': xreopen(optarg, "w", stdout); break;
            default: return 1;
        }
    }
    
    if(optind + 3 > argc){
        fprintf(stderr, "Usage: bwsplice samse [-n max_occ] [-f out.sam] [-r RG_line] <prefix> <in.sai> <in.fq>\n");
        return 1;
    }
    bwt_sai2sam_se_core(argv[optind], argv[optind+1], argv[optind + 2], n_occ);
    free(bwt_rg_line);
    free(bwt_rg_id);
    return 0;
}
