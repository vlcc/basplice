#ifndef BWASE_H
#define BWASE_H

#include "2BWT-Interface.h"
#include "bwtaln.h"

#ifdef __cplusplus
extern "C" {
#endif

	// Initialize mapping tables in the bwa single-end mapper.
	void bwase_initialize();
	// Calculate the approximate position of the sequence from the specified bwt with loaded suffix array.
	void bwa_cal_pac_pos_core(const Idx2BWT *bi_bwt, bwa_seq_t* seq, const int max_mm, const float fnr);
	// Refine the approximate position of the sequence to an actual placement for the sequence.
	void bwa_refine_gapped(const HSP *hsp, int n_seqs, bwa_seq_t *seqs, HSP *hsp_cs);
	// Backfill certain alignment properties mainly centering around number of matches.
	void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s);
    void bwt_aln2seq_splicing(const Idx2BWT *bi_bwt, bwa_seq_t *seq, int n_aln, bwt_aln1_t *aln);
	// Calculate the end position of a read given a certain sequence.
	int64_t pos_end(const bwa_seq_t *p);
	//
	//bwtint_t bwa_sa2pos(const bntseq_t *bns, const bwt_t *bwt, bwtint_t sapos, int len, int *strand);

#ifdef __cplusplus
}
#endif

#endif // BWASE_H
