#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

gap_opt_t *gap_init_opt()
{
	gap_opt_t *o;
	o = (gap_opt_t*)calloc(2, sizeof(gap_opt_t));
	/* IMPORTANT: s_mm*10 should be about the average base error
	   rate. Voilating this requirement will break pairing! */
	o->s_mm = 3;
    o->s_gapo = 11;
    o->s_gape = 4;
	o->max_diff = -1;
    o->max_gapo = 1;
    o->max_gape = 6;
	o->indel_end_skip = 5;
    o->max_del_occ = 10;
    o->max_entries = 2000000;
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	o->seed_len = 32;
    o->max_seed_diff = 2;
	o->fnr = 0.04;
	o->n_threads = 1;
	o->max_top2 = 30;
	o->trim_qual = 0;
	return o;
}

int bwa_cal_maxdiff(int l, double err, double thres)
{
	double elambda = exp(-l * err);
	double sum, y = 1.0;
	int k, x = 1;
	for (k = 1, sum = elambda; k < 1000; ++k) {
		y *= l * err;
		x *= k;
		sum += elambda * y / x;
		if (1.0 - sum < thres) return k;
	}
	return 2;
}

void get_seq_from_ref(const Idx2BWT *bi_bwt, bwtint_t _pos, int len, ubyte_t *ref_seq)
{
    bwtint_t k;
    HSP *hsp = bi_bwt->hsp;
    int l = 0;
    for(k = _pos; k < _pos + len; k ++)
        ref_seq[l++] = hsp->packedDNA[k>>4] >>((~k & 15) <<  1) & 3;
    return ref_seq;
}

// width must be filled as zero
// str is origin chars, not converted
// 0 : foreward, 1 - backward
int bwt_cal_width(const Idx2BWT *bi_bwt, int len, const ubyte_t *str, bwt_width_t *width,
        int type)
{
	unsigned int k, l;
    //ubyte_t *conv_str;
	int i, bid;
	bid = 0;
	k = 0;
    l = bi_bwt->bwt->textLength;

    //convert str to 2BWT recognised coding scheme
    //BWTConvertPattern(bi_bwt, str, len, str);
    //printf("%s\n", str);

    if(type == 0){
	    for (i = len - 1; i > 0; --i) {
	    	ubyte_t c = str[i];
	    	if (c < 4) {
                BWTSARangeBackward(bi_bwt, c, &k, &l);
	    	}
	    	if (k > l || c > 3) { // then restart
	    		k = 0;
	    		l = bi_bwt->bwt->textLength;
	    		++bid;
	    	}
	    	width[i].w = l - k + 1;
	    	width[i].bid = bid;
	    }
    } else {
        for(i = 0 ; i < len; ++i){
            ubyte_t c = str[i];
            if(c < 4)
                BWTSARangeForeward(bi_bwt, c, &k , &l);
            if(k > l || c >3){
                k = 0;
                l = bi_bwt->bwt->textLength;
                ++bid;
            }
            width[i].w = l -k + 1;
            width[i].bid = bid;
        }
    }
	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}
// for debug, get seq real position from seq name
void get_pos_from_name(char const *name, bwtint_t *start, bwtint_t *end)
{
    size_t pos, len;
    char *tmp_c, *token, *tmp;
    int i, swp;
    len = strlen(name);
    tmp_c = calloc(len + 1, sizeof(char));
    memset(tmp_c, 0, len+1);
    memcpy(tmp_c, name, len);
    tmp = tmp_c;
    for(token = strtok(tmp,"_"), i = 0;
            token != NULL;
            token = strtok(NULL, "_"), ++i){
        if(i == 1){
            *start = atoi(token);
        }
        if(i == 2) {
            *end = atoi(token);
            break;
        }
    }
    if (*start > *end) {
        swp = *start;
        *start = *end;
        *end = swp;
    }
    free(tmp_c);
}

/* // if seq_id & offset is found, return 1, else insert into pos
 * int inline bwt_check_pos(pos_t *pos, int seq_id, unsigned int offset, int n_seq_id){
 *     int i, j, k, m, n;
 *     offset_t *q, *p;
 *     for (n = 0; n < n_seq_id && pos->seq_id != 0; ++n){
 *         if(pos->seq_id == seq_id){
 *             q = pos->offset;
 *             while(q != NULL){
 *                 if(q->offset > offset){
 *                     //the list is ordered
 *                     p = (offset_t *)malloc(sizeof(offset));
 *                     p->next = q->next;
 *                     p->offset = offset;
 *                     q->next = p;
 *                     return 0;
 *                 }else if(q->offset == offset){
 *                     return 1;
 *                 }else{
 *                     q = q->next;
 *                 }
 *             }
 *         }
 *         else
 *             pos += 1;
 *     }
 *     pos->seq_id = seq_id;
 *     p = (offset_t *)malloc(sizeof(offset));
 *     p->offset = offset;
 *     pos->offset = p;
 *     return 0
 * }
 * void bwt_seed_aln_correlative_check(bwt_seed_aln_t *multi_seed_aln,
 *         int seed_num)
 * {
 *     int i, j, n, n_seq_id, seq_id;
 *     bwt_aln1_t *aln;
 *     bwtint_t k, l, rev_k, rev_l, m;
 *     unsigned int offset;
 *     bwt_seed_aln_t *p;
 *     pos_t *pos;
 *     n_seq_id = 30;
 *     pos = (pos_t *)calloc(n_seq_id, sizeof(pos_t));
 * 
 *     for(i = 0; i < seed_num; ++i){
 *         p = multi_seed_aln + i;
 *         for(j = 0; j < p->n_aln; ++j){
 *             aln = p + j;
 *             for(m = 0, m < aln->l - aln->l + 1; ++m){
 *                 //get seq_ids and offsets for an alignment, by saIndex k - l
 *                 BWTRetrievePositionFromSAIndex(bi_bwt, aln->k + m, &seq_id, &offset);
 *             }
 *         }
 *     }
 * }
 */

// 去除多余的aln，找出没有的aln上的部分
//TODO ust bit operation as another implement 
/* bwt_list_t *bwt_extend_result_check2(bwt_seed_aln_t *multi_seed_aln,
 *         int seed_num, int len)
 * {
 *     int i, *mapped_info,j;
 *     bwt_aln1_t *tmp_aln;
 *     bwt_list_t *q, *head, *x;
 *     mapped_info = (int *)calloc(len + 1, sizeof(int));
 *     for(i=0; i < seed_num; ++i){
 *         tmp_aln = (multi_seed_aln+i)->aln;
 *         for(j = tmp_aln->start ; j <= tmp_aln->end; ++j)
 *             mapped_info[j]=1;
 *     }
 *     i = 1;
 *     head = NULL;
 *     while(i <= len){
 *         if(mapped_info[i] == 0){
 *             j = i;
 *             while(mapped_info[j]==0 && j <= len)
 *                 ++j;
 *             x = (bwt_list_t *)malloc(sizeof(bwt_list_t));
 *             x->start = i;
 *             x->end = j - 1;
 *             x->next = NULL;
 *             if(head == NULL){
 *                 q = x;
 *                 head = x;
 *             }else{
 *                 q->next = x;
 *                 q = q->next;
 *             }
 *             i = j;
 *         }else{
 *             while(mapped_info[i] == 1 && i <= len)
 *                 ++i;
 *         }
 *     }
 *     return head;
 * }
 */

//seqs是所有需要mapping的reads
void bwa_cal_sa_reg_gap(int tid,const Idx2BWT *bi_bwt, int n_seqs, bwa_seq_t *seqs,
        const gap_opt_t *opt, bwt_array_t *arr)
{
	int i, max_l = 0, max_len;
    //char tmp_c;
    //unsigned int offset;
	gap_stack_t *stack;
	bwt_width_t *width_seed, *width_back, *width_fore;
	gap_opt_t local_opt = *opt;
    bwt_aux_t *aux;

    // init aux struct for mapping
    aux = (bwt_aux_t *)malloc(sizeof(bwt_aux_t));
    aux->bi_bwt = bi_bwt;
    aux->opt = opt;

    // bwt array
    aux->arr = arr;

    //get the max_len of all seqs
	for (i = max_len = 0; i != n_seqs; ++i)
		if (seqs[i].len > max_len)
            max_len = seqs[i].len;
    aux->max_len = max_len;

    //get max_diff base on error rate 
	if (opt->fnr > 0.0)
        local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo)
        local_opt.max_gapo = local_opt.max_diff;

    // init stack for bracktracing
	aux->stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);


    // init bwt_width_t 
    aux->width_back = calloc(max_len+1, sizeof(bwt_width_t));
    aux->width_fore = calloc(max_len+1, sizeof(bwt_width_t));
    aux->width_seed = calloc(max_len+1, sizeof(bwt_width_t));

    // init reverse complimetn seq
    aux->rc_seq = (ubyte_t *)calloc(max_len + 1, sizeof(ubyte_t));

    // TODO for debug
    bwtint_t seq_id, offset, hit_num, seq_start, seq_end, called_num;
    bwtint_t miss_num, tmp_sa;
    int b_hit = 0;
    called_num = 0; hit_num = 0; offset = 0; seq_id = 0; seq_id = 0;
    offset = 0; seq_start = 0; seq_end = 0, miss_num = 0;

    int _i, j, _j;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
#ifdef HAVE_PTHREAD
		if (i % opt->n_threads != tid)
            continue;
#endif

        // check whether too many N in seq
        for(j = _j = 0; j < p->len; ++j)
            if(p->seq[j] > 3) ++ _j;
        if(_j > local_opt.max_diff)
            continue;

		p->sa = 0;
        p->type = BWA_TYPE_NO_MATCH;
        p->c1 = p->c2 = 0;
        p->n_aln = 0;
        p->aln = 0;
        aux->seq = p->seq;
        aux->len = p->len;

        // cal max diff and set seed length
		if (opt->fnr > 0.0)
            aux->opt->max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
		aux->opt->seed_len = opt->seed_len < p->len? opt->seed_len : 0x7fffffff;

        // init rc_seq
        memset(aux->rc_seq, 0, max_len * sizeof(ubyte_t));
        memcpy(aux->rc_seq, p->seq, p->len * sizeof(ubyte_t));
        seq_reverse(p->len, aux->rc_seq, 1);

        // reset width_back
	    memset(aux->width_back, 0, (max_len + 1) * sizeof(bwt_width_t));

        // reverse complement seq first
        for(_i = 1; _i >= 0; --_i) {
            if(p->len > aux->opt->seed_len)
		    	bwt_cal_width(bi_bwt, opt->seed_len, (_i == 0? p->seq : aux->rc_seq)+
                        (p->len - aux->opt->seed_len) , aux->width_seed,1);

		    bwt_cal_width(bi_bwt, p->len,_i == 0? p->seq: aux->rc_seq, aux->width_back, 1);
            aux->strand = _i;
		    p->aln = bwt_match_gap(aux, &p->n_aln);
            if(p->n_aln != 0) {
                for(_j = 0; _j < p->n_aln; ++_j)
                    (p->aln+_j)->strand = _i;
                break;
            }
        }

        // 判断是否mapping上，如没有则做splicing mapping
        if(p->n_aln == 0){
            aux->opt = &local_opt;
            p->aln = bwt_splice_match(aux, &p->n_aln);
        } else{
            p->aln->start = 0;   //能mapping到ref上的reads起始位置
            p->aln->end = p->len - 1;
        }
        // TODO debug
        if(p->n_aln != 0) {
            b_hit = 0;
            called_num ++;
            get_pos_from_name(p->name, &seq_start, &seq_end);
            for(tmp_sa = p->aln->k; tmp_sa <= p->aln->l; ++tmp_sa){
                BWTRetrievePositionFromSAIndex(bi_bwt, tmp_sa, &seq_id, &offset);
                if(seq_start == 0 && seq_end ==0)
                    continue;
                if(offset >= seq_start && offset < seq_end) {
                    hit_num ++;
                    b_hit = 1;
                    break;
                }
            }
            if(b_hit == 0) {
                miss_num += 1;
            }
        }

		// clean up the unused data in the record
		if(p->name) free(p->name);
        if(p->seq)  free(p->seq);
        if(p->rseq) free(p->rseq);
        if(p->qual) free(p->qual);
		p->name = 0;
        p->seq = p->rseq = p->qual = 0;

	}
	free(aux->width_seed);
	free(aux->width_fore);
    free(aux->width_back);
    free(aux->rc_seq);
	gap_destroy_stack(aux->stack);
    free(aux);
    fprintf(stderr, "call seq num: %u, hit seq num: %u\n", called_num, hit_num);
}

#ifdef HAVE_PTHREAD
typedef struct {
	int tid;
	Idx2BWT *bi_bwt;
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
    const bwt_array_t *arr;
} thread_aux_t;

static void *worker(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	bwa_cal_sa_reg_gap(d->tid, d->bi_bwt, d->n_seqs, d->seqs, d->opt, d->arr);
	return 0;
}
#endif

bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa)
{
	bwa_seqio_t *ks;
	if (mode & BWA_MODE_BAM) { // open BAM
		int which = 0;
		if (mode & BWA_MODE_BAM_SE) which |= 4;
		if (mode & BWA_MODE_BAM_READ1) which |= 1;
		if (mode & BWA_MODE_BAM_READ2) which |= 2;
		if (which == 0) which = 7; // then read all reads
		ks = bwa_bam_open(fn_fa, which);
	} else ks = bwa_seq_open(fn_fa);
	return ks;
}

void bwa_aln_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt)
{
	int i, n_seqs, tot_seqs = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t;
	Idx2BWT *bi_bwt;

	// initialization
	ks = bwa_open_reads(opt->mode, fn_fa);

	{ // load 2BWT
		char *str = (char*)calloc(strlen(prefix) + 10, 1);
		strcpy(str, prefix); 
        strcat(str, ".index");  
        bi_bwt = BWTLoad2BWT(str,".sa");
		free(str);
	}

    // init array to record splicing site
    bwt_array_t *arr;
    arr = bwt_array_init();

	// core loop
	err_fwrite(opt, sizeof(gap_opt_t), 1, stdout);
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
		tot_seqs += n_seqs;
		t = clock();
		fprintf(stderr, "[bwa_aln_core] calculate SA coordinate... \n");

#ifdef HAVE_PTHREAD
		if (opt->n_threads <= 1) { // no multi-threading at all
			bwa_cal_sa_reg_gap(0, bi_bwt, n_seqs, seqs, opt, arr);
		} else {
			pthread_t *tid;
			pthread_attr_t attr;
			thread_aux_t *data;
			int j;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
			tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
			for (j = 0; j < opt->n_threads; ++j) {
				data[j].tid = j; data[j].bi_bwt = bi_bwt;
				data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
                data[j].arr = arr;
				pthread_create(&tid[j], &attr, worker, data + j);
			}
			for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
			free(data); free(tid);
		}
#else

		bwa_cal_sa_reg_gap(0, bi_bwt, n_seqs, seqs, opt, arr);
#endif


		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
		t = clock();
		fprintf(stderr, "[bwa_aln_core] write to the disk... ");
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
/*             fprintf(stderr,"n_aln: %d,%d\n", p->n_aln, i);
 */
			err_fwrite(&p->n_aln, sizeof(int), 1, stdout);
			if (p->n_aln) err_fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, stdout);
		}
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		bwa_free_read_seq(n_seqs, seqs);
		fprintf(stderr, "[bwa_aln_core] %d sequences have been processed.\n", tot_seqs);
	}

	// destroy
	BWTFree2BWT(bi_bwt);
	bwa_seq_close(ks);
    bwt_array_destory(arr);
}

int bwa_aln(int argc, char *argv[])
{
	int c, opte = -1;
	gap_opt_t *opt;

	opt = gap_init_opt();
	while ((c = getopt(argc, argv, "n:o:e:i:d:l:k:cLR:m:t:NM:O:E:q:f:b012IYB:")) >= 0) {
		switch (c) {
		case 'n':
			if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
			else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
			break;
		case 'o': opt->max_gapo = atoi(optarg); break;
		case 'e': opte = atoi(optarg); break;
		case 'M': opt->s_mm = atoi(optarg); break;
		case 'O': opt->s_gapo = atoi(optarg); break;
		case 'E': opt->s_gape = atoi(optarg); break;
		case 'd': opt->max_del_occ = atoi(optarg); break;
		case 'i': opt->indel_end_skip = atoi(optarg); break;
		case 'l': opt->seed_len = atoi(optarg); break;
		case 'k': opt->max_seed_diff = atoi(optarg); break;
		case 'm': opt->max_entries = atoi(optarg); break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
		case 'R': opt->max_top2 = atoi(optarg); break;
		case 'q': opt->trim_qual = atoi(optarg); break;
		case 'c': opt->mode &= ~BWA_MODE_COMPREAD; break;
		case 'N': opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff; break;
		case 'f': xreopen(optarg, "wb", stdout); break;
		case 'b': opt->mode |= BWA_MODE_BAM; break;
		case '0': opt->mode |= BWA_MODE_BAM_SE; break;
		case '1': opt->mode |= BWA_MODE_BAM_READ1; break;
		case '2': opt->mode |= BWA_MODE_BAM_READ2; break;
		case 'I': opt->mode |= BWA_MODE_IL13; break;
		case 'Y': opt->mode |= BWA_MODE_CFY; break;
		case 'B': opt->mode |= atoi(optarg) << 24; break;
		default: return 1;
		}
	}
	if (opte > 0) {
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa aln [options] <prefix> <in.fq>\n\n");
		fprintf(stderr, "Options: -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n",
				BWA_AVG_ERR, opt->fnr);
		fprintf(stderr, "         -o INT    maximum number or fraction of gap opens [%d]\n", opt->max_gapo);
		fprintf(stderr, "         -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
		fprintf(stderr, "         -i INT    do not put an indel within INT bp towards the ends [%d]\n", opt->indel_end_skip);
		fprintf(stderr, "         -d INT    maximum occurrences for extending a long deletion [%d]\n", opt->max_del_occ);
		fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
		fprintf(stderr, "         -k INT    maximum differences in the seed [%d]\n", opt->max_seed_diff);
		fprintf(stderr, "         -m INT    maximum entries in the queue [%d]\n", opt->max_entries);
		fprintf(stderr, "         -t INT    number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "         -M INT    mismatch penalty [%d]\n", opt->s_mm);
		fprintf(stderr, "         -O INT    gap open penalty [%d]\n", opt->s_gapo);
		fprintf(stderr, "         -E INT    gap extension penalty [%d]\n", opt->s_gape);
		fprintf(stderr, "         -R INT    stop searching when there are >INT equally best hits [%d]\n", opt->max_top2);
		fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
        fprintf(stderr, "         -f FILE   file to write output to instead of stdout\n");
		fprintf(stderr, "         -B INT    length of barcode\n");
//		fprintf(stderr, "         -c        input sequences are in the color space\n");
		fprintf(stderr, "         -L        log-scaled gap penalty for long deletions\n");
		fprintf(stderr, "         -N        non-iterative mode: search for all n-difference hits (slooow)\n");
		fprintf(stderr, "         -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
		fprintf(stderr, "         -b        the input read file is in the BAM format\n");
		fprintf(stderr, "         -0        use single-end reads only (effective with -b)\n");
		fprintf(stderr, "         -1        use the 1st read in a pair (effective with -b)\n");
		fprintf(stderr, "         -2        use the 2nd read in a pair (effective with -b)\n");
		fprintf(stderr, "         -Y        filter Casava-filtered sequences\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (opt->fnr > 0.0) {
		int i, k;
		for (i = 17, k = 0; i <= 250; ++i) {
			int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, opt->fnr);
			if (l != k) fprintf(stderr, "[bwa_aln] %dbp reads: max_diff = %d\n", i, l);
			k = l;
		}
	}
	bwa_aln_core(argv[optind], argv[optind+1], opt);
	free(opt);
	return 0;
}

/* rgoya: Temporary clone of aln_path2cigar to accomodate for bwa_cigar_t,
__cigar_op and __cigar_len while keeping stdaln stand alone */
bwa_cigar_t *bwa_aln_path2cigar(const path_t *path, int path_len, int *n_cigar)
{
	uint32_t *cigar32;
	bwa_cigar_t *cigar;
	int i;
	cigar32 = aln_path2cigar32((path_t*) path, path_len, n_cigar);
	cigar = (bwa_cigar_t*)cigar32;
	for (i = 0; i < *n_cigar; ++i)
                cigar[i] = __cigar_create( (cigar32[i]&0xf), (cigar32[i]>>4) );
	return cigar;
}
