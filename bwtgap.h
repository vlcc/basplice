#ifndef __BWTGAP_H_
#define __BWTGAP_H_

#include "2BWT-Interface.h"
#include "bwtaln.h"

typedef struct list_t{
    int start;
    int end;
    struct list_t *next;
}bwt_list_t;

typedef struct{
    int n_aln;
    int strand;
    bwt_aln1_t *aln;
}bwt_seed_aln_t;

#ifdef __cplusplus
extern "C" {
#endif

    void seq_reverse(int len, ubyte_t *seq, int is_comp);
	gap_stack_t *gap_init_stack(int max_mm, int max_gapo, int max_gape, const gap_opt_t *opt);
	void gap_destroy_stack(gap_stack_t *stack);
    bwt_aln1_t *bwt_match_gap(bwt_aux_t *aux, int *_n_aln);
	void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s);

    int bwt_extend_backward(bwt_aux_t *aux, bwt_aln1_t *aln, int *_left);
    int bwt_extend_foreward(bwt_aux_t *aux, bwt_aln1_t *aln, int *_right);
    bwt_aln1_t *bwt_splice_match(bwt_aux_t *aux, int *_n_aln);

#ifdef __cplusplus
}
#endif

#endif
