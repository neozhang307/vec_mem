#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include "kstring.h"
#include "bwamem.h"
#include "bntseq.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

/* Theory on probability and scoring *ungapped* alignment
 *
 * s'(a,b) = log[P(b|a)/P(b)] = log[4P(b|a)], assuming uniform base distribution
 * s'(a,a) = log(4), s'(a,b) = log(4e/3), where e is the error rate
 *
 * Scale s'(a,b) to s(a,a) s.t. s(a,a)=x. Then s(a,b) = x*s'(a,b)/log(4), or conversely: s'(a,b)=s(a,b)*log(4)/x
 *
 * If the matching score is x and mismatch penalty is -y, we can compute error rate e:
 *   e = .75 * exp[-log(4) * y/x]
 *
 * log P(seq) = \sum_i log P(b_i|a_i) = \sum_i {s'(a,b) - log(4)}
 *   = \sum_i { s(a,b)*log(4)/x - log(4) } = log(4) * (S/x - l)
 *
 * where S=\sum_i s(a,b) is the alignment score. Converting to the phred scale:
 *   Q(seq) = -10/log(10) * log P(seq) = 10*log(4)/log(10) * (l - S/x) = 6.02 * (l - S/x)
 *
 *
 * Gap open (zero gap): q' = log[P(gap-open)], r' = log[P(gap-ext)] (see Durbin et al. (1998) Section 4.1)
 * Then q = x*log[P(gap-open)]/log(4), r = x*log[P(gap-ext)]/log(4)
 *
 * When there are gaps, l should be the length of alignment matches (i.e. the M operator in CIGAR)
 */

static const bntseq_t *global_bns = 0; // for debugging only

mem_opt_t *mem_opt_init()
{
	mem_opt_t *o;
	o = calloc(1, sizeof(mem_opt_t));
	o->flag = 0;
	o->a = 1; o->b = 4;
	o->o_del = o->o_ins = 6;
	o->e_del = o->e_ins = 1;
	o->w = 100;
	o->T = 30;
	o->zdrop = 100;
	o->pen_unpaired = 17;
	o->pen_clip5 = o->pen_clip3 = 5;

	o->max_mem_intv = 20;

	o->min_seed_len = 19;
	o->split_width = 10;
	o->max_occ = 500;
	o->max_chain_gap = 10000;
	o->max_ins = 10000;
	o->mask_level = 0.50;
	o->drop_ratio = 0.50;
	o->XA_drop_ratio = 0.80;
	o->split_factor = 1.5;
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->max_XA_hits = 5;
	o->max_XA_hits_alt = 200;
	o->max_matesw = 50;
	o->mask_level_redun = 0.95;
	o->min_chain_weight = 0;
	o->max_chain_extend = 1<<30;
	o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
	bwa_fill_scmat(o->a, o->b, o->mat);
	return o;
}

/***************************
 * Collection SA invervals *
 ***************************/

#define intv_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mem_intv, bwtintv_t, intv_lt)

typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

static smem_aux_t *smem_aux_init()
{
	smem_aux_t *a;
	a = calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = calloc(1, sizeof(bwtintv_v));
	return a;
}

static void smem_aux_destroy(smem_aux_t *a)
{	
	free(a->tmpv[0]->a); free(a->tmpv[0]);
	free(a->tmpv[1]->a); free(a->tmpv[1]);
	free(a->mem.a); free(a->mem1.a);
	free(a);
}

static void mem_collect_intv(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a)
{
	int i, k, x = 0, old_n;
	int start_width = 1;
	int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
	a->mem.n = 0;
	// first pass: find all SMEMs
	while (x < len) {
		if (seq[x] < 4) {
			x = bwt_smem1(bwt, len, seq, x, start_width, &a->mem1, a->tmpv);
			for (i = 0; i < a->mem1.n; ++i) {
				bwtintv_t *p = &a->mem1.a[i];
				int slen = (uint32_t)p->info - (p->info>>32); // seed length
				if (slen >= opt->min_seed_len)
					kv_push(bwtintv_t, a->mem, *p);
			}
		} else ++x;
	}
	// second pass: find MEMs inside a long SMEM
	old_n = a->mem.n;
	for (k = 0; k < old_n; ++k) {
		bwtintv_t *p = &a->mem.a[k];
		int start = p->info>>32, end = (int32_t)p->info;
		if (end - start < split_len || p->x[2] > opt->split_width) continue;
		bwt_smem1(bwt, len, seq, (start + end)>>1, p->x[2]+1, &a->mem1, a->tmpv);
		for (i = 0; i < a->mem1.n; ++i)
			if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info>>32) >= opt->min_seed_len)
				kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
	}
	// third pass: LAST-like
	if (opt->max_mem_intv > 0) {
		x = 0;
		while (x < len) {
			if (seq[x] < 4) {
				if (1) {
					bwtintv_t m;
					x = bwt_seed_strategy1(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
					if (m.x[2] > 0) kv_push(bwtintv_t, a->mem, m);
				} else { // for now, we never come to this block which is slower
					x = bwt_smem1a(bwt, len, seq, x, start_width, opt->max_mem_intv, &a->mem1, a->tmpv);
					for (i = 0; i < a->mem1.n; ++i)
						kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
				}
			} else ++x;
		}
	}
	// sort
	ks_introsort(mem_intv, a->mem.n, a->mem.a);
}

/************
 * Chaining *
 ************/

typedef struct {
	int64_t rbeg;
	int32_t qbeg, len;
	int score;
} mem_seed_t; // unaligned memory

typedef struct {
	int n, m, first, rid;
	uint32_t w:29, kept:2, is_alt:1;
	float frac_rep;
	int64_t pos;
	mem_seed_t *seeds;
} mem_chain_t;

typedef struct { size_t n, m; mem_chain_t *a;  } mem_chain_v;

#include "kbtree.h"

#define chain_cmp(a, b) (((b).pos < (a).pos) - ((a).pos < (b).pos))
KBTREE_INIT(chn, mem_chain_t, chain_cmp)

// return 1 if the seed is merged into the chain
static int test_and_merge(const mem_opt_t *opt, int64_t l_pac, mem_chain_t *c, const mem_seed_t *p, int seed_rid)
{
	int64_t qend, rend, x, y;
	const mem_seed_t *last = &c->seeds[c->n-1];
	qend = last->qbeg + last->len;
	rend = last->rbeg + last->len;
	if (seed_rid != c->rid) return 0; // different chr; request a new chain
	if (p->qbeg >= c->seeds[0].qbeg && p->qbeg + p->len <= qend && p->rbeg >= c->seeds[0].rbeg && p->rbeg + p->len <= rend)
		return 1; // contained seed; do nothing
	if ((last->rbeg < l_pac || c->seeds[0].rbeg < l_pac) && p->rbeg >= l_pac) return 0; // don't chain if on different strand
	x = p->qbeg - last->qbeg; // always non-negtive
	y = p->rbeg - last->rbeg;
	if (y >= 0 && x - y <= opt->w && y - x <= opt->w && x - last->len < opt->max_chain_gap && y - last->len < opt->max_chain_gap) { // grow the chain
		if (c->n == c->m) {
			c->m <<= 1;
			c->seeds = realloc(c->seeds, c->m * sizeof(mem_seed_t));
		}
		c->seeds[c->n++] = *p;
		return 1;
	}
	return 0; // request to add a new chain
}

int mem_chain_weight(const mem_chain_t *c)
{
	int64_t end;
	int j, w = 0, tmp;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->qbeg >= end) w += s->len;
		else if (s->qbeg + s->len > end) w += s->qbeg + s->len - end;
		end = end > s->qbeg + s->len? end : s->qbeg + s->len;
	}
	tmp = w; w = 0;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->rbeg >= end) w += s->len;
		else if (s->rbeg + s->len > end) w += s->rbeg + s->len - end;
		end = end > s->rbeg + s->len? end : s->rbeg + s->len;
	}
	w = w < tmp? w : tmp;
	return w < 1<<30? w : (1<<30)-1;
}

void mem_print_chain(const bntseq_t *bns, mem_chain_v *chn)
{
	int i, j;
	for (i = 0; i < chn->n; ++i) {
		mem_chain_t *p = &chn->a[i];
		err_printf("* Found CHAIN(%d): n=%d; weight=%d", i, p->n, mem_chain_weight(p));
		for (j = 0; j < p->n; ++j) {
			bwtint_t pos;
			int is_rev;
			pos = bns_depos(bns, p->seeds[j].rbeg, &is_rev);
			if (is_rev) pos -= p->seeds[j].len - 1;
			err_printf("\t%d;%d;%d,%ld(%s:%c%ld)", p->seeds[j].score, p->seeds[j].len, p->seeds[j].qbeg, (long)p->seeds[j].rbeg, bns->anns[p->rid].name, "+-"[is_rev], (long)(pos - bns->anns[p->rid].offset) + 1);
		}
		err_putchar('\n');
	}
}

mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, void *buf)
{
	int i, b, e, l_rep;
	int64_t l_pac = bns->l_pac;
	mem_chain_v chain;
	kbtree_t(chn) *tree;
	smem_aux_t *aux;

	kv_init(chain);
	if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match
	tree = kb_init(chn, KB_DEFAULT_SIZE);

	aux = buf? (smem_aux_t*)buf : smem_aux_init();
	mem_collect_intv(opt, bwt, len, seq, aux);
	for (i = 0, b = e = l_rep = 0; i < aux->mem.n; ++i) { // compute frac_rep
		bwtintv_t *p = &aux->mem.a[i];
		int sb = (p->info>>32), se = (uint32_t)p->info;
		if (p->x[2] <= opt->max_occ) continue;
		if (sb > e) l_rep += e - b, b = sb, e = se;
		else e = e > se? e : se;
	}
	l_rep += e - b;
	for (i = 0; i < aux->mem.n; ++i) {
		bwtintv_t *p = &aux->mem.a[i];
		int step, count, slen = (uint32_t)p->info - (p->info>>32); // seed length
		int64_t k;
		// if (slen < opt->min_seed_len) continue; // ignore if too short or too repetitive
		step = p->x[2] > opt->max_occ? p->x[2] / opt->max_occ : 1;
		for (k = count = 0; k < p->x[2] && count < opt->max_occ; k += step, ++count) {
			mem_chain_t tmp, *lower, *upper;
			mem_seed_t s;
			int rid, to_add = 0;
			s.rbeg = tmp.pos = bwt_sa(bwt, p->x[0] + k); // this is the base coordinate in the forward-reverse reference
			s.qbeg = p->info>>32;
			s.score= s.len = slen;
			rid = bns_intv2rid(bns, s.rbeg, s.rbeg + s.len);
			if (rid < 0) continue; // bridging multiple reference sequences or the forward-reverse boundary; TODO: split the seed; don't discard it!!!
			if (kb_size(tree)) {
				kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain
				if (!lower || !test_and_merge(opt, l_pac, lower, &s, rid)) to_add = 1;
			} else to_add = 1;
			if (to_add) { // add the seed as a new chain
				tmp.n = 1; tmp.m = 4;
				tmp.seeds = calloc(tmp.m, sizeof(mem_seed_t));
				tmp.seeds[0] = s;
				tmp.rid = rid;
				tmp.is_alt = !!bns->anns[rid].is_alt;
				kb_putp(chn, tree, &tmp);
			}
		}
	}
	if (buf == 0) smem_aux_destroy(aux);

	kv_resize(mem_chain_t, chain, kb_size(tree));

	#define traverse_func(p_) (chain.a[chain.n++] = *(p_))
	__kb_traverse(mem_chain_t, tree, traverse_func);
	#undef traverse_func

	for (i = 0; i < chain.n; ++i) chain.a[i].frac_rep = (float)l_rep / len;
	if (bwa_verbose >= 4) printf("* fraction of repetitive seeds: %.3f\n", (float)l_rep / len);

	kb_destroy(chn, tree);
	return chain;
}

/********************
 * Filtering chains *
 ********************/

#define chn_beg(ch) ((ch).seeds->qbeg)
#define chn_end(ch) ((ch).seeds[(ch).n-1].qbeg + (ch).seeds[(ch).n-1].len)

#define flt_lt(a, b) ((a).w > (b).w)
KSORT_INIT(mem_flt, mem_chain_t, flt_lt)

int mem_chain_flt(const mem_opt_t *opt, int n_chn, mem_chain_t *a)
{
	int i, k;
	kvec_t(int) chains = {0,0,0}; // this keeps int indices of the non-overlapping chains
	if (n_chn == 0) return 0; // no need to filter
	// compute the weight of each chain and drop chains with small weight
	for (i = k = 0; i < n_chn; ++i) {
		mem_chain_t *c = &a[i];
		c->first = -1; c->kept = 0;
		c->w = mem_chain_weight(c);
		if (c->w < opt->min_chain_weight) free(c->seeds);
		else a[k++] = *c;
	}
	n_chn = k;
	ks_introsort(mem_flt, n_chn, a);
	// pairwise chain comparisons
	a[0].kept = 3;
	kv_push(int, chains, 0);
	for (i = 1; i < n_chn; ++i) {
		int large_ovlp = 0;
		for (k = 0; k < chains.n; ++k) {
			int j = chains.a[k];
			int b_max = chn_beg(a[j]) > chn_beg(a[i])? chn_beg(a[j]) : chn_beg(a[i]);
			int e_min = chn_end(a[j]) < chn_end(a[i])? chn_end(a[j]) : chn_end(a[i]);
			if (e_min > b_max && (!a[j].is_alt || a[i].is_alt)) { // have overlap; don't consider ovlp where the kept chain is ALT while the current chain is primary
				int li = chn_end(a[i]) - chn_beg(a[i]);
				int lj = chn_end(a[j]) - chn_beg(a[j]);
				int min_l = li < lj? li : lj;
				if (e_min - b_max >= min_l * opt->mask_level && min_l < opt->max_chain_gap) { // significant overlap
					large_ovlp = 1;
					if (a[j].first < 0) a[j].first = i; // keep the first shadowed hit s.t. mapq can be more accurate
					if (a[i].w < a[j].w * opt->drop_ratio && a[j].w - a[i].w >= opt->min_seed_len<<1)
						break;
				}
			}
		}
		if (k == chains.n) {
			kv_push(int, chains, i);
			a[i].kept = large_ovlp? 2 : 3;
		}
	}
	for (i = 0; i < chains.n; ++i) {
		mem_chain_t *c = &a[chains.a[i]];
		if (c->first >= 0) a[c->first].kept = 1;
	}
	free(chains.a);
	for (i = k = 0; i < n_chn; ++i) { // don't extend more than opt->max_chain_extend .kept=1/2 chains
		if (a[i].kept == 0 || a[i].kept == 3) continue;
		if (++k >= opt->max_chain_extend) break;
	}
	for (; i < n_chn; ++i)
		if (a[i].kept < 3) a[i].kept = 0;
	for (i = k = 0; i < n_chn; ++i) { // free discarded chains
		mem_chain_t *c = &a[i];
		if (c->kept == 0) free(c->seeds);
		else a[k++] = a[i];
	}
	return k;
}

/******************************
 * De-overlap single-end hits *
 ******************************/

#define alnreg_slt2(a, b) ((a).re < (b).re)
KSORT_INIT(mem_ars2, mem_alnreg_t, alnreg_slt2)

#define alnreg_slt(a, b) ((a).score > (b).score || ((a).score == (b).score && ((a).rb < (b).rb || ((a).rb == (b).rb && (a).qb < (b).qb))))
KSORT_INIT(mem_ars, mem_alnreg_t, alnreg_slt)

#define alnreg_hlt(a, b)  ((a).score > (b).score || ((a).score == (b).score && ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash, mem_alnreg_t, alnreg_hlt)

#define alnreg_hlt2(a, b) ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && ((a).score > (b).score || ((a).score == (b).score && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash2, mem_alnreg_t, alnreg_hlt2)

#define PATCH_MAX_R_BW 0.05f
#define PATCH_MIN_SC_RATIO 0.90f

int mem_patch_reg(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, const mem_alnreg_t *a, const mem_alnreg_t *b, int *_w)
{
	int w, score, q_s, r_s;
	double r;
	if (bns == 0 || pac == 0 || query == 0) return 0;
	assert(a->rid == b->rid && a->rb <= b->rb);
	if (a->rb < bns->l_pac && b->rb >= bns->l_pac) return 0; // on different strands
	if (a->qb >= b->qb || a->qe >= b->qe || a->re >= b->re) return 0; // not colinear
	w = (a->re - b->rb) - (a->qe - b->qb); // required bandwidth
	w = w > 0? w : -w; // l = abs(l)
	r = (double)(a->re - b->rb) / (b->re - a->rb) - (double)(a->qe - b->qb) / (b->qe - a->qb); // relative bandwidth
	r = r > 0.? r : -r; // r = fabs(r)
	if (bwa_verbose >= 4)
		printf("* potential hit merge between [%d,%d)<=>[%ld,%ld) and [%d,%d)<=>[%ld,%ld), @ %s; w=%d, r=%.4g\n",
			   a->qb, a->qe, (long)a->rb, (long)a->re, b->qb, b->qe, (long)b->rb, (long)b->re, bns->anns[a->rid].name, w, r);
	if (a->re < b->rb || a->qe < b->qb) { // no overlap on query or on ref
		if (w > opt->w<<1 || r >= PATCH_MAX_R_BW) return 0; // the bandwidth or the relative bandwidth is too large
	} else if (w > opt->w<<2 || r >= PATCH_MAX_R_BW*2) return 0; // more permissive if overlapping on both ref and query
	// global alignment
	w += a->w + b->w;
	w = w < opt->w<<2? w : opt->w<<2;
	if (bwa_verbose >= 4) printf("* test potential hit merge with global alignment; w=%d\n", w);
	bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w, bns->l_pac, pac, b->qe - a->qb, query + a->qb, a->rb, b->re, &score, 0, 0);
	q_s = (int)((double)(b->qe - a->qb) / ((b->qe - b->qb) + (a->qe - a->qb)) * (b->score + a->score) + .499); // predicted score from query
	r_s = (int)((double)(b->re - a->rb) / ((b->re - b->rb) + (a->re - a->rb)) * (b->score + a->score) + .499); // predicted score from ref
	if (bwa_verbose >= 4) printf("* score=%d;(%d,%d)\n", score, q_s, r_s);
	if ((double)score / (q_s > r_s? q_s : r_s) < PATCH_MIN_SC_RATIO) return 0;
	*_w = w;
	return score;
}

int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, int n, mem_alnreg_t *a)
{
	int m, i, j;
	if (n <= 1) return n;
	ks_introsort(mem_ars2, n, a); // sort by the END position, not START!
	for (i = 0; i < n; ++i) a[i].n_comp = 1;
	for (i = 1; i < n; ++i) {
		mem_alnreg_t *p = &a[i];
		if (p->rid != a[i-1].rid || p->rb >= a[i-1].re + opt->max_chain_gap) continue; // then no need to go into the loop below
		for (j = i - 1; j >= 0 && p->rid == a[j].rid && p->rb < a[j].re + opt->max_chain_gap; --j) {
			mem_alnreg_t *q = &a[j];
			int64_t or, oq, mr, mq;
			int score, w;
			if (q->qe == q->qb) continue; // a[j] has been excluded
			or = q->re - p->rb; // overlap length on the reference
			oq = q->qb < p->qb? q->qe - p->qb : p->qe - q->qb; // overlap length on the query
			mr = q->re - q->rb < p->re - p->rb? q->re - q->rb : p->re - p->rb; // min ref len in alignment
			mq = q->qe - q->qb < p->qe - p->qb? q->qe - q->qb : p->qe - p->qb; // min qry len in alignment
			if (or > opt->mask_level_redun * mr && oq > opt->mask_level_redun * mq) { // one of the hits is redundant
				if (p->score < q->score) {
					p->qe = p->qb;
					break;
				} else q->qe = q->qb;
			} else if (q->rb < p->rb && (score = mem_patch_reg(opt, bns, pac, query, q, p, &w)) > 0) { // then merge q into p
				p->n_comp += q->n_comp + 1;
				p->seedcov = p->seedcov > q->seedcov? p->seedcov : q->seedcov;
				p->sub = p->sub > q->sub? p->sub : q->sub;
				p->csub = p->csub > q->csub? p->csub : q->csub;
				p->qb = q->qb, p->rb = q->rb;
				p->truesc = p->score = score;
				p->w = w;
				q->qb = q->qe;
			}
		}
	}
	for (i = 0, m = 0; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	n = m;
	ks_introsort(mem_ars, n, a);
	for (i = 1; i < n; ++i) { // mark identical hits
		if (a[i].score == a[i-1].score && a[i].rb == a[i-1].rb && a[i].qb == a[i-1].qb)
			a[i].qe = a[i].qb;
	}
	for (i = 1, m = 1; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	return m;
}

typedef kvec_t(int) int_v;

static void mem_mark_primary_se_core(const mem_opt_t *opt, int n, mem_alnreg_t *a, int_v *z)
{ // similar to the loop in mem_chain_flt()
	int i, k, tmp;
	tmp = opt->a + opt->b;
	tmp = opt->o_del + opt->e_del > tmp? opt->o_del + opt->e_del : tmp;
	tmp = opt->o_ins + opt->e_ins > tmp? opt->o_ins + opt->e_ins : tmp;
	z->n = 0;
	kv_push(int, *z, 0);
	for (i = 1; i < n; ++i) {
		for (k = 0; k < z->n; ++k) {
			int j = z->a[k];
			int b_max = a[j].qb > a[i].qb? a[j].qb : a[i].qb;
			int e_min = a[j].qe < a[i].qe? a[j].qe : a[i].qe;
			if (e_min > b_max) { // have overlap
				int min_l = a[i].qe - a[i].qb < a[j].qe - a[j].qb? a[i].qe - a[i].qb : a[j].qe - a[j].qb;
				if (e_min - b_max >= min_l * opt->mask_level) { // significant overlap
					if (a[j].sub == 0) a[j].sub = a[i].score;
					if (a[j].score - a[i].score <= tmp && (a[j].is_alt || !a[i].is_alt))
						++a[j].sub_n;
					break;
				}
			}
		}
		if (k == z->n) kv_push(int, *z, i);
		else a[i].secondary = z->a[k];
	}
}

int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id)
{
	int i, n_pri;
	int_v z = {0,0,0};
	if (n == 0) return 0;
	for (i = n_pri = 0; i < n; ++i) {
		a[i].sub = a[i].alt_sc = 0, a[i].secondary = a[i].secondary_all = -1, a[i].hash = hash_64(id+i);
		if (!a[i].is_alt) ++n_pri;
	}
	ks_introsort(mem_ars_hash, n, a);
	mem_mark_primary_se_core(opt, n, a, &z);
	for (i = 0; i < n; ++i) {
		mem_alnreg_t *p = &a[i];
		p->secondary_all = i; // keep the rank in the first round
		if (!p->is_alt && p->secondary >= 0 && a[p->secondary].is_alt)
			p->alt_sc = a[p->secondary].score;
	}
	if (n_pri >= 0 && n_pri < n) {
		kv_resize(int, z, n);
		if (n_pri > 0) ks_introsort(mem_ars_hash2, n, a);
		for (i = 0; i < n; ++i) z.a[a[i].secondary_all] = i;
		for (i = 0; i < n; ++i) {
			if (a[i].secondary >= 0) {
				a[i].secondary_all = z.a[a[i].secondary];
				if (a[i].is_alt) a[i].secondary = INT_MAX;
			} else a[i].secondary_all = -1;
		}
		if (n_pri > 0) { // mark primary for hits to the primary assembly only
			for (i = 0; i < n_pri; ++i) a[i].sub = 0, a[i].secondary = -1;
			mem_mark_primary_se_core(opt, n_pri, a, &z);
		}
	} else {
		for (i = 0; i < n; ++i)
			a[i].secondary_all = a[i].secondary;
	}
	free(z.a);
	return n_pri;
}

/*********************************
 * Test if a seed is good enough *
 *********************************/

#define MEM_SHORT_EXT 50
#define MEM_SHORT_LEN 200

#define MEM_HSP_COEF 1.1f
#define MEM_MINSC_COEF 5.5f
#define MEM_SEEDSW_COEF 0.05f

int mem_seed_sw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_seed_t *s)
{
	int qb, qe, rid;
	int64_t rb, re, mid, l_pac = bns->l_pac;
	uint8_t *rseq = 0;
	kswr_t x;

	if (s->len >= MEM_SHORT_LEN) return -1; // the seed is longer than the max-extend; no need to do SW
	qb = s->qbeg, qe = s->qbeg + s->len;
	rb = s->rbeg, re = s->rbeg + s->len;
	mid = (rb + re) >> 1;
	qb -= MEM_SHORT_EXT; qb = qb > 0? qb : 0;
	qe += MEM_SHORT_EXT; qe = qe < l_query? qe : l_query;
	rb -= MEM_SHORT_EXT; rb = rb > 0? rb : 0;
	re += MEM_SHORT_EXT; re = re < l_pac<<1? re : l_pac<<1;
	if (rb < l_pac && l_pac < re) {
		if (mid < l_pac) re = l_pac;
		else rb = l_pac;
	}
	if (qe - qb >= MEM_SHORT_LEN || re - rb >= MEM_SHORT_LEN) return -1; // the seed seems good enough; no need to do SW

	rseq = bns_fetch_seq(bns, pac, &rb, mid, &re, &rid);
	x = ksw_align2(qe - qb, (uint8_t*)query + qb, re - rb, rseq, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, KSW_XSTART, 0);
	free(rseq);
	return x.score;
}

void mem_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, int n_chn, mem_chain_t *a)
{
	double min_l = opt->min_chain_weight? MEM_HSP_COEF * opt->min_chain_weight : MEM_MINSC_COEF * log(l_query);
	int i, j, k, min_HSP_score = (int)(opt->a * min_l + .499);
	if (min_l > MEM_SEEDSW_COEF * l_query) return; // don't run the following for short reads
	for (i = 0; i < n_chn; ++i) {
		mem_chain_t *c = &a[i];
		for (j = k = 0; j < c->n; ++j) {
			mem_seed_t *s = &c->seeds[j];
			s->score = mem_seed_sw(opt, bns, pac, l_query, query, s);
			if (s->score < 0 || s->score >= min_HSP_score) {
				s->score = s->score < 0? s->len * opt->a : s->score;
				c->seeds[k++] = *s;
			}
		}
		c->n = k;
	}
}

/****************************************
 * Construct the alignment from a chain *
 ****************************************/
//NEO:
//max mismatch no larger than opt->w<<1
//a is the match score
static inline int cal_max_gap(const mem_opt_t *opt, int qlen)
{
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}

#define MAX_BAND_TRY  2
//#define DEBUG

static long count = 0;
void add_count(long data)
{
    count+=data;
}
void reset()
{
    count=0;
}
void printcount()
{
    fprintf(stderr,"the count is %ld now\n",count);
}

void mem_chain2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av)
{
	int i, k, rid, max_off[2], aw[2]; // aw: actual bandwidth used in extension
	int64_t l_pac = bns->l_pac, rmax[2], tmp, max = 0;
	const mem_seed_t *s;
	uint8_t *rseq = 0;
	uint64_t *srt;

    
    /*
     NEO:
     @para rmax[2] {thread private}:
                [0]: begin
                [1]: end
     */
    
	if (c->n == 0) return;
	// get the max possible span
    // NEO: should set on CPU
    // NEO: can we orgnize query by max possible span?
	rmax[0] = l_pac<<1; rmax[1] = 0;
	for (i = 0; i < c->n; ++i) {
		int64_t b, e;
		const mem_seed_t *t = &c->seeds[i];
		b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
		e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
		rmax[0] = rmax[0] < b? rmax[0] : b;
		rmax[1] = rmax[1] > e? rmax[1] : e;
		if (t->len > max) max = t->len;
	}
	rmax[0] = rmax[0] > 0? rmax[0] : 0;
	rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
	if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
		if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
		else rmax[0] = l_pac;
	}
#ifdef DEBUG
    add_count(rmax[1]-rmax[0]);
#endif
	// retrieve the reference sequence
	rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);//NEO: potentially OOM, in every 10MB batch, average 67MB
	assert(c->rid == rid);
    
    // NEO:
    // external sorting
    // Generate str: str is an index,   high 32 bit is SW score (sorted)
    //                                  low 32 bit is real index
	srt = malloc(c->n * 8);
	for (i = 0; i < c->n; ++i)
		srt[i] = (uint64_t)c->seeds[i].score<<32 | i;
	ks_introsort_64(c->n, srt);// NEO: srt in decending order

    
    // NEO: should do modification in this part in the future
	for (k = c->n - 1; k >= 0; --k) {
		mem_alnreg_t *a;
		s = &c->seeds[(uint32_t)srt[k]];

        
        // NEO: this part is belong to CPU, should migrate this to the end of this function.
        // NEO:
        // should know how many seed would be drop in this place
        // Test if the seed is in future align
		for (i = 0; i < av->n; ++i) { // test whether extension has been made before
			mem_alnreg_t *p = &av->a[i];
			int64_t rd;
			int qd, w, max_gap;
			if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue; // not fully contained
			if (s->len - p->seedlen0 > .1 * l_query) continue; // this seed may give a better alignment
			// qd: distance ahead of the seed on query; rd: on reference
			qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
			max_gap = cal_max_gap(opt, qd < rd? qd : rd); // the maximal gap allowed in regions ahead of the seed
			w = max_gap < p->w? max_gap : p->w; // bounded by the band width
			if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
			// similar to the previous four lines, but this time we look at the region behind
			qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
			max_gap = cal_max_gap(opt, qd < rd? qd : rd);
			w = max_gap < p->w? max_gap : p->w;
			if (qd - rd < w && rd - qd < w) break;
		}
        // NEO:
        // rescue the seed marked as overlap, if it would lead to a different result
		if (i < av->n) { // the seed is (almost) contained in an existing alignment; further testing is needed to confirm it is not leading to a different aln
			if (bwa_verbose >= 4)
				printf("** Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment [%d,%d) <=> [%ld,%ld)\n",
					   k, (long)s->len, (long)s->qbeg, (long)s->rbeg, av->a[i].qb, av->a[i].qe, (long)av->a[i].rb, (long)av->a[i].re);
            
            //NEO: block structure
			for (i = k + 1; i < c->n; ++i) { // check overlapping seeds in the same chain
				const mem_seed_t *t;
				if (srt[i] == 0) continue;
				t = &c->seeds[(uint32_t)srt[i]];
				if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODO: more efficient by early stopping
				if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
				if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
			}
			
            
            if (i == c->n) { // no overlapping seeds; then skip extension
				srt[k] = 0; // mark that seed extension has not been performed
				continue;
			}
			if (bwa_verbose >= 4)
				printf("** Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", k);
		}
        
		a = kv_pushp(mem_alnreg_t, *av);
		memset(a, 0, sizeof(mem_alnreg_t));
		a->w = aw[0] = aw[1] = opt->w;
		a->score = a->truesc = -1;
		a->rid = c->rid;

		if (bwa_verbose >= 4) err_printf("** ---> Extending from seed(%d) [%ld;%ld,%ld] @ %s <---\n", k, (long)s->len, (long)s->qbeg, (long)s->rbeg, bns->anns[c->rid].name);
		if (s->qbeg) { // left extension
			uint8_t *rs, *qs;
			int qle, tle, gtle, gscore;
			qs = malloc(s->qbeg);
			for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
			tmp = s->rbeg - rmax[0];
			rs = malloc(tmp);
			for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[0] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Left ref:   "); for (j = 0; j < tmp; ++j) putchar("ACGTN"[(int)rs[j]]); putchar('\n');
					printf("*** Left query: "); for (j = 0; j < s->qbeg; ++j) putchar("ACGTN"[(int)qs[j]]); putchar('\n');
				}
                //NEO: the most time consuming part
				a->score = ksw_extend2(s->qbeg, qs, tmp, rs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, s->len * opt->a, &qle, &tle, &gtle, &gscore, &max_off[0]);
				if (bwa_verbose >= 4) { printf("*** Left extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[0], max_off[0]); fflush(stdout); }
				if (a->score == prev || max_off[0] < (aw[0]>>1) + (aw[0]>>2)) break;
			}
			// check whether we prefer to reach the end of the query
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
				a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
				a->truesc = a->score;
			} else { // to-end extension
				a->qb = 0, a->rb = s->rbeg - gtle;
				a->truesc = gscore;
			}
			free(qs); free(rs);
		} else a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;

		if (s->qbeg + s->len != l_query) { // right extension
			int qle, tle, qe, re, gtle, gscore, sc0 = a->score;
			qe = s->qbeg + s->len;
			re = s->rbeg + s->len - rmax[0];
			assert(re >= 0);
            
            //NEO: warp or block
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[1] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Right ref:   "); for (j = 0; j < rmax[1] - rmax[0] - re; ++j) putchar("ACGTN"[(int)rseq[re+j]]); putchar('\n');
					printf("*** Right query: "); for (j = 0; j < l_query - qe; ++j) putchar("ACGTN"[(int)query[qe+j]]); putchar('\n');
				}
				a->score = ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
				if (bwa_verbose >= 4) { printf("*** Right extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[1], max_off[1]); fflush(stdout); }
				if (a->score == prev || max_off[1] < (aw[1]>>1) + (aw[1]>>2)) break;
			}
            
			// similar to the above
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
				a->qe = qe + qle, a->re = rmax[0] + re + tle;
				a->truesc += a->score - sc0;
			} else { // to-end extension
				a->qe = l_query, a->re = rmax[0] + re + gtle;
				a->truesc += gscore - sc0;
			}
		} else a->qe = l_query, a->re = s->rbeg + s->len;
		if (bwa_verbose >= 4) printf("*** Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}_bandwidth={%d,%d}\n", a->qb, a->qe, (long)a->rb, (long)a->re, a->score, aw[0], aw[1]);

		// compute seedcov
		for (i = 0, a->seedcov = 0; i < c->n; ++i) {
			const mem_seed_t *t = &c->seeds[i];
			if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
				a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
		}
		a->w = aw[0] > aw[1]? aw[0] : aw[1];
		a->seedlen0 = s->len;

		a->frac_rep = c->frac_rep;
	}
	free(srt); free(rseq);
}


/*****************************
 * Basic hit->SAM conversion *
 *****************************/

static inline int infer_bw(int l1, int l2, int score, int a, int q, int r)
{
	int w;
	if (l1 == l2 && l1 * a - score < (q + r - a)<<1) return 0; // to get equal alignment length, we need at least two gaps
	w = ((double)((l1 < l2? l1 : l2) * a - score - q) / r + 2.);
	if (w < abs(l1 - l2)) w = abs(l1 - l2);
	return w;
}

static inline int get_rlen(int n_cigar, const uint32_t *cigar)
{
	int k, l;
	for (k = l = 0; k < n_cigar; ++k) {
		int op = cigar[k]&0xf;
		if (op == 0 || op == 2)
			l += cigar[k]>>4;
	}
	return l;
}

void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_)
{
	int i, l_name;
	mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert

	if (m_) mtmp = *m_, m = &mtmp;
	// set flag
	p->flag |= m? 0x1 : 0; // is paired in sequencing
	p->flag |= p->rid < 0? 0x4 : 0; // is mapped
	p->flag |= m && m->rid < 0? 0x8 : 0; // is mate mapped
	if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
		p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
	if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
		m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
	p->flag |= p->is_rev? 0x10 : 0; // is on the reverse strand
	p->flag |= m && m->is_rev? 0x20 : 0; // is mate on the reverse strand

	// print up to CIGAR
	l_name = strlen(s->name);
	ks_resize(str, str->l + s->l_seq + l_name + (s->qual? s->l_seq : 0) + 20);
	kputsn(s->name, l_name, str); kputc('\t', str); // QNAME
	kputw((p->flag&0xffff) | (p->flag&0x10000? 0x100 : 0), str); kputc('\t', str); // FLAG
	if (p->rid >= 0) { // with coordinate
		kputs(bns->anns[p->rid].name, str); kputc('\t', str); // RNAME
		kputl(p->pos + 1, str); kputc('\t', str); // POS
		kputw(p->mapq, str); kputc('\t', str); // MAPQ
		if (p->n_cigar) { // aligned
			for (i = 0; i < p->n_cigar; ++i) {
				int c = p->cigar[i]&0xf;
				if (!(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt && (c == 3 || c == 4))
					c = which? 4 : 3; // use hard clipping for supplementary alignments
				kputw(p->cigar[i]>>4, str); kputc("MIDSH"[c], str);
			}
		} else kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
	} else kputsn("*\t0\t0\t*", 7, str); // without coordinte
	kputc('\t', str);

	// print the mate position if applicable
	if (m && m->rid >= 0) {
		if (p->rid == m->rid) kputc('=', str);
		else kputs(bns->anns[m->rid].name, str);
		kputc('\t', str);
		kputl(m->pos + 1, str); kputc('\t', str);
		if (p->rid == m->rid) {
			int64_t p0 = p->pos + (p->is_rev? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
			int64_t p1 = m->pos + (m->is_rev? get_rlen(m->n_cigar, m->cigar) - 1 : 0);
			if (m->n_cigar == 0 || p->n_cigar == 0) kputc('0', str);
			else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
		} else kputc('0', str);
	} else kputsn("*\t0\t0", 5, str);
	kputc('\t', str);

	// print SEQ and QUAL
	if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
		kputsn("*\t*", 3, str);
	} else if (!p->is_rev) { // the forward strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) { // have cigar && not the primary alignment && not softclip all
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qb += p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qe -= p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qb; i < qe; ++i) str->s[str->l++] = "ACGTN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qb; i < qe; ++i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	} else { // the reverse strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) {
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qe -= p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qb += p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qe-1; i >= qb; --i) str->s[str->l++] = "TGCAN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qe-1; i >= qb; --i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	}

	// print optional tags
	if (p->n_cigar) {
		kputsn("\tNM:i:", 6, str); kputw(p->NM, str);
		kputsn("\tMD:Z:", 6, str); kputs((char*)(p->cigar + p->n_cigar), str);
	}
	if (p->score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p->score, str); }
	if (p->sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(p->sub, str); }
	if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
	if (!(p->flag & 0x100)) { // not multi-hit
		for (i = 0; i < n; ++i)
			if (i != which && !(list[i].flag&0x100)) break;
		if (i < n) { // there are other primary hits; output them
			kputsn("\tSA:Z:", 6, str);
			for (i = 0; i < n; ++i) {
				const mem_aln_t *r = &list[i];
				int k;
				if (i == which || (r->flag&0x100)) continue; // proceed if: 1) different from the current; 2) not shadowed multi hit
				kputs(bns->anns[r->rid].name, str); kputc(',', str);
				kputl(r->pos+1, str); kputc(',', str);
				kputc("+-"[r->is_rev], str); kputc(',', str);
				for (k = 0; k < r->n_cigar; ++k) {
					kputw(r->cigar[k]>>4, str); kputc("MIDSH"[r->cigar[k]&0xf], str);
				}
				kputc(',', str); kputw(r->mapq, str);
				kputc(',', str); kputw(r->NM, str);
				kputc(';', str);
			}
		}
		if (p->alt_sc > 0)
			ksprintf(str, "\tpa:f:%.3f", (double)p->score / p->alt_sc);
	}
	if (p->XA) { kputsn("\tXA:Z:", 6, str); kputs(p->XA, str); }
	if (s->comment) { kputc('\t', str); kputs(s->comment, str); }
	if ((opt->flag&MEM_F_REF_HDR) && p->rid >= 0 && bns->anns[p->rid].anno != 0 && bns->anns[p->rid].anno[0] != 0) {
		int tmp;
		kputsn("\tXR:Z:", 6, str);
		tmp = str->l;
		kputs(bns->anns[p->rid].anno, str);
		for (i = tmp; i < str->l; ++i) // replace TAB in the comment to SPACE
			if (str->s[i] == '\t') str->s[i] = ' ';
	}
	kputc('\n', str);
}

/************************
 * Integrated interface *
 ************************/

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a)
{
	int mapq, l, sub = a->sub? a->sub : opt->min_seed_len * opt->a;
	double identity;
	sub = a->csub > sub? a->csub : sub;
	if (sub >= a->score) return 0;
	l = a->qe - a->qb > a->re - a->rb? a->qe - a->qb : a->re - a->rb;
	identity = 1. - (double)(l * opt->a - a->score) / (opt->a + opt->b) / l;
	if (a->score == 0) {
		mapq = 0;
	} else if (opt->mapQ_coef_len > 0) {
		double tmp;
		tmp = l < opt->mapQ_coef_len? 1. : opt->mapQ_coef_fac / log(l);
		tmp *= identity * identity;
		mapq = (int)(6.02 * (a->score - sub) / opt->a * tmp * tmp + .499);
	} else {
		mapq = (int)(MEM_MAPQ_COEF * (1. - (double)sub / a->score) * log(a->seedcov) + .499);
		mapq = identity < 0.95? (int)(mapq * identity * identity + .499) : mapq;
	}
	if (a->sub_n > 0) mapq -= (int)(4.343 * log(a->sub_n+1) + .499);
	if (mapq > 60) mapq = 60;
	if (mapq < 0) mapq = 0;
	mapq = (int)(mapq * (1. - a->frac_rep) + .499);
	return mapq;
}

// TODO (future plan): group hits into a uint64_t[] array. This will be cleaner and more flexible
void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m)
{
	extern char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, mem_alnreg_v *a, int l_query, const char *query);
	kstring_t str;
	kvec_t(mem_aln_t) aa;
	int k, l;
	char **XA = 0;

	if (!(opt->flag & MEM_F_ALL))
		XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq);
	kv_init(aa);
	str.l = str.m = 0; str.s = 0;
	for (k = l = 0; k < a->n; ++k) {
		mem_alnreg_t *p = &a->a[k];
		mem_aln_t *q;
		if (p->score < opt->T) continue;
		if (p->secondary >= 0 && (p->is_alt || !(opt->flag&MEM_F_ALL))) continue;
		if (p->secondary >= 0 && p->secondary < INT_MAX && p->score < a->a[p->secondary].score * opt->drop_ratio) continue;
		q = kv_pushp(mem_aln_t, aa);
		*q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p);
		assert(q->rid >= 0); // this should not happen with the new code
		q->XA = XA? XA[k] : 0;
		q->flag |= extra_flag; // flag secondary
		if (p->secondary >= 0) q->sub = -1; // don't output sub-optimal score
		if (l && p->secondary < 0) // if supplementary
			q->flag |= (opt->flag&MEM_F_NO_MULTI)? 0x10000 : 0x800;
		if (l && !p->is_alt && q->mapq > aa.a[0].mapq) q->mapq = aa.a[0].mapq;
		++l;
	}
	if (aa.n == 0) { // no alignments good enough; then write an unaligned record
		mem_aln_t t;
		t = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, 0);
		t.flag |= extra_flag;
		mem_aln2sam(opt, bns, &str, s, 1, &t, 0, m);
	} else {
		for (k = 0; k < aa.n; ++k)
			mem_aln2sam(opt, bns, &str, s, aa.n, aa.a, k, m);
		for (k = 0; k < aa.n; ++k) free(aa.a[k].cigar);
		free(aa.a);
	}
	s->sam = str.s;
	if (XA) {
		for (k = 0; k < a->n; ++k) free(XA[k]);
		free(XA);
	}
}


mem_alnreg_v mem_align1_core(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, void *buf)
{
	int i;
	mem_chain_v chn;
	mem_alnreg_v regs;

	for (i = 0; i < l_seq; ++i) // convert to 2-bit encoding if we have not done so
		seq[i] = seq[i] < 4? seq[i] : nst_nt4_table[(int)seq[i]];

	chn = mem_chain(opt, bwt, bns, l_seq, (uint8_t*)seq, buf);
	chn.n = mem_chain_flt(opt, chn.n, chn.a);
	mem_flt_chained_seeds(opt, bns, pac, l_seq, (uint8_t*)seq, chn.n, chn.a);
	if (bwa_verbose >= 4) mem_print_chain(bns, &chn);

	kv_init(regs);
	for (i = 0; i < chn.n; ++i) {
		mem_chain_t *p = &chn.a[i];
		if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", i);
		mem_chain2aln(opt, bns, pac, l_seq, (uint8_t*)seq, p, &regs);
		free(chn.a[i].seeds);
	}
	free(chn.a);
 //   fprintf(stderr, "count chain_v: %ld\n",g_count1);
 //   fprintf(stderr, "count ndrop: %ld\n",g_count2);
	regs.n = mem_sort_dedup_patch(opt, bns, pac, (uint8_t*)seq, regs.n, regs.a);
	if (bwa_verbose >= 4) {
		err_printf("* %ld chains remain after removing duplicated chains\n", regs.n);
		for (i = 0; i < regs.n; ++i) {
			mem_alnreg_t *p = &regs.a[i];
			printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
		}
	}
	for (i = 0; i < regs.n; ++i) {
		mem_alnreg_t *p = &regs.a[i];
		if (p->rid >= 0 && bns->anns[p->rid].is_alt)
			p->is_alt = 1;
	}
	return regs;
}

mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const char *query_, const mem_alnreg_t *ar)
{
	mem_aln_t a;
	int i, w2, tmp, qb, qe, NM, score, is_rev, last_sc = -(1<<30), l_MD;
	int64_t pos, rb, re;
	uint8_t *query;

	memset(&a, 0, sizeof(mem_aln_t));
	if (ar == 0 || ar->rb < 0 || ar->re < 0) { // generate an unmapped record
		a.rid = -1; a.pos = -1; a.flag |= 0x4;
		return a;
	}
	qb = ar->qb, qe = ar->qe;
	rb = ar->rb, re = ar->re;
	query = malloc(l_query);
	for (i = 0; i < l_query; ++i) // convert to the nt4 encoding
		query[i] = query_[i] < 5? query_[i] : nst_nt4_table[(int)query_[i]];
	a.mapq = ar->secondary < 0? mem_approx_mapq_se(opt, ar) : 0;
	if (ar->secondary >= 0) a.flag |= 0x100; // secondary alignment
	tmp = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_del, opt->e_del);
	w2  = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_ins, opt->e_ins);
	w2 = w2 > tmp? w2 : tmp;
	if (bwa_verbose >= 4) printf("* Band width: inferred=%d, cmd_opt=%d, alnreg=%d\n", w2, opt->w, ar->w);
	if (w2 > opt->w) w2 = w2 < ar->w? w2 : ar->w;
	i = 0; a.cigar = 0;
	do {
		free(a.cigar);
		w2 = w2 < opt->w<<2? w2 : opt->w<<2;
		a.cigar = bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w2, bns->l_pac, pac, qe - qb, (uint8_t*)&query[qb], rb, re, &score, &a.n_cigar, &NM);
		if (bwa_verbose >= 4) printf("* Final alignment: w2=%d, global_sc=%d, local_sc=%d\n", w2, score, ar->truesc);
		if (score == last_sc || w2 == opt->w<<2) break; // it is possible that global alignment and local alignment give different scores
		last_sc = score;
		w2 <<= 1;
	} while (++i < 3 && score < ar->truesc - opt->a);
	l_MD = strlen((char*)(a.cigar + a.n_cigar)) + 1;
	a.NM = NM;
	pos = bns_depos(bns, rb < bns->l_pac? rb : re - 1, &is_rev);
	a.is_rev = is_rev;
	if (a.n_cigar > 0) { // squeeze out leading or trailing deletions
		if ((a.cigar[0]&0xf) == 2) {
			pos += a.cigar[0]>>4;
			--a.n_cigar;
			memmove(a.cigar, a.cigar + 1, a.n_cigar * 4 + l_MD);
		} else if ((a.cigar[a.n_cigar-1]&0xf) == 2) {
			--a.n_cigar;
			memmove(a.cigar + a.n_cigar, a.cigar + a.n_cigar + 1, l_MD); // MD needs to be moved accordingly
		}
	}
	if (qb != 0 || qe != l_query) { // add clipping to CIGAR
		int clip5, clip3;
		clip5 = is_rev? l_query - qe : qb;
		clip3 = is_rev? qb : l_query - qe;
		a.cigar = realloc(a.cigar, 4 * (a.n_cigar + 2) + l_MD);
		if (clip5) {
			memmove(a.cigar+1, a.cigar, a.n_cigar * 4 + l_MD); // make room for 5'-end clipping
			a.cigar[0] = clip5<<4 | 3;
			++a.n_cigar;
		}
		if (clip3) {
			memmove(a.cigar + a.n_cigar + 1, a.cigar + a.n_cigar, l_MD); // make room for 3'-end clipping
			a.cigar[a.n_cigar++] = clip3<<4 | 3;
		}
	}
	a.rid = bns_pos2rid(bns, pos);
	assert(a.rid == ar->rid);
	a.pos = pos - bns->anns[a.rid].offset;
	a.score = ar->score; a.sub = ar->sub > ar->csub? ar->sub : ar->csub;
	a.is_alt = ar->is_alt; a.alt_sc = ar->alt_sc;
	free(query);
	return a;
}

typedef struct {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	const mem_pestat_t *pes;
	smem_aux_t **aux;
	bseq1_t *seqs;
	mem_alnreg_v *regs;
	int64_t n_processed;
} worker_t;



static void worker1(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		if (bwa_verbose >= 4) printf("=====> Processing read '%s' <=====\n", w->seqs[i].name);
		w->regs[i] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[i].l_seq, w->seqs[i].seq, w->aux[tid]);
	} else {
		if (bwa_verbose >= 4) printf("=====> Processing read '%s'/1 <=====\n", w->seqs[i<<1|0].name);
		w->regs[i<<1|0] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[i<<1|0].l_seq, w->seqs[i<<1|0].seq, w->aux[tid]);
		if (bwa_verbose >= 4) printf("=====> Processing read '%s'/2 <=====\n", w->seqs[i<<1|1].name);
		w->regs[i<<1|1] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[i<<1|1].l_seq, w->seqs[i<<1|1].seq, w->aux[tid]);
	}
}

static void worker2(void *data, int i, int tid)
{
	extern int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]);
	extern void mem_reg2ovlp(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a);
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		if (bwa_verbose >= 4) printf("=====> Finalizing read '%s' <=====\n", w->seqs[i].name);
		mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a, w->n_processed + i);
		mem_reg2sam(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
		free(w->regs[i].a);
	} else {
		if (bwa_verbose >= 4) printf("=====> Finalizing read pair '%s' <=====\n", w->seqs[i<<1|0].name);
		mem_sam_pe(w->opt, w->bns, w->pac, w->pes, (w->n_processed>>1) + i, &w->seqs[i<<1], &w->regs[i<<1]);
		free(w->regs[i<<1|0].a); free(w->regs[i<<1|1].a);
	}
}

/*********************************************************/
/*this segment is added by Lingqi Zhang*/
void mem_chain2maxspan(const mem_opt_t *opt, const mem_chain_t *c, int l_query, int64_t l_pac,  int64_t rmax[2])
{
    /*
     NEO:
     @para rmax[2] {thread private}:
     [0]: begin
     [1]: end
     */
    // get the max possible span
    // NEO: should set on CPU
    // NEO: can we orgnize query by max possible span?
    // NEO: potentially OOM, in every 10MB batch, average 67MB
    int i;
    int64_t max=0;
    rmax[0] = l_pac<<1; rmax[1] = 0;
    for (i = 0; i < c->n; ++i) {
        int64_t b, e;
        const mem_seed_t *t = &c->seeds[i];
        b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
        e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
        rmax[0] = rmax[0] < b? rmax[0] : b;
        rmax[1] = rmax[1] > e? rmax[1] : e;
        if (t->len > max) max = t->len;
    }
    rmax[0] = rmax[0] > 0? rmax[0] : 0;
    rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
    if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
        if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
        else rmax[0] = l_pac;
    }
#ifdef DEBUG
    add_count(rmax[1]-rmax[0]);
#endif
}
#include"ksw_batch_simd.h"
//typedef struct
//{
//    int qlen;
//    const uint8_t *query;
//    int rlen;
//    const uint8_t *ref;
//}swseq_t;
//typedef struct
//{
//    int score;
//    int max_off;
//    
//    int qle,tle;
//    int gtle, gscore;
//    
//    int h0;
//    swseq_t* sw_seq;
//}swrst_t;

void mem_chain2aln_preextent(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const uint8_t *query_rev, int64_t rmax[2], const uint8_t *rseq, const uint8_t *rseq_rev, const mem_chain_t *c,  swseq_t * forward, swseq_t * backward, swrst_t* swfwd_, swrst_t* swbwd_)
{
    int k;
    int64_t  tmp;
    const mem_seed_t *s;
    
    for (k = c->n - 1; k >= 0; --k) {
        s = &c->seeds[k];
        swrst_t *swfwd = &swfwd_[k];
        swrst_t *swbwd = &swbwd_[k];
        swseq_t *swf = &forward[k];
        swseq_t *swb = &backward[k];
        memset(swf, 0, sizeof(swseq_t));
        memset(swb, 0, sizeof(swseq_t));

        swfwd->sw_seq = swf;
        swbwd->sw_seq = swb;
        swfwd->score =  s->len * opt->a;
        if (s->qbeg) { // left extension
            const uint8_t  *qs;
            const uint8_t *rs;

            qs = &query_rev[l_query-s->qbeg];
            
            tmp = s->rbeg - rmax[0];
            rs = &rseq_rev[rmax[1]-s->rbeg];
            swf->qlen = s->qbeg;
            swf->query = qs;
            swf->rlen = tmp;
            swf->ref = rs;
            
//            swfwd->score =s->len * opt->a;
            swfwd->h0 =s->len * opt->a;
        }
        swbwd->score =  s->len * opt->a;
        if (s->qbeg + s->len != l_query) { // right extension
            int qe, re;// sc0 = swfwd->score;
            qe = s->qbeg + s->len;
            re = s->rbeg + s->len - rmax[0];
            assert(re >= 0);
            
            swb->qlen = l_query - qe;
            swb->query = query + qe;
            swb->rlen = rmax[1] - rmax[0] - re;
            swb->ref = rseq + re;
            //swbwd->score = sc0;
        }
    }
}

void mem_chain_extent(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, swrst_t* sws)
{
    const mem_seed_t *s;
    
    for (int k = 0; k < c->n; ++k) {
        s = &c->seeds[k];
        swrst_t *sw = &sws[k];
        swseq_t *seq = sw->sw_seq;
        if(seq->qlen!=0)
        {
            sw->score = ksw_extend2_mod(seq->qlen, seq->query, seq->rlen,seq->ref, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, opt->zdrop, sw->h0, &sw->qle, &sw->tle, &sw->gtle, &sw->gscore, &sw->max_off);
        }
        else{
            sw->score =  s->len * opt->a;
        }
    }
}


void mem_chain2aln_swextent(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, swrst_t* swfwd_, swrst_t* swbwd_)
{
    int k;
    
    //left extend
    mem_chain_extent(opt, bns, pac, l_query,query, c, swfwd_);
    //init left extend
    for (k = 0; k < c->n; ++k) {
        int pre_score= swfwd_[k].score;
        swrst_t *swbwd = &swbwd_[k];
        swbwd->h0=pre_score;
    }
    //right extend
    mem_chain_extent(opt, bns, pac, l_query,query, c, swbwd_);
}

void mem_chain2aln_post(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, int64_t rmax[2], mem_alnreg_t *av_firstpass, swrst_t* swfwd_, swrst_t* swbwd_)
{
    int i, k;//, aw; // aw: actual bandwidth used in extension
    const mem_seed_t *s;
    
    
    for (k = c->n - 1; k >= 0; --k) {
        s = &c->seeds[k];
        
        mem_alnreg_t *a;
        a = &av_firstpass[k];
        memset(a, 0, sizeof(mem_alnreg_t));
      //  a->w = aw = 0;//opt->w;
        a->score = a->truesc = -1;
        a->rid = c->rid;
        swrst_t *swfwd = &swfwd_[k];
        swrst_t *swbwd = &swbwd_[k];
        swseq_t *fwd = swfwd->sw_seq;
        swseq_t *bwd = swbwd->sw_seq;
        if(fwd->qlen!=0)//left extension
        {
            int qle, tle, gtle, gscore;
            
            a->score = swfwd->score;
            qle = swfwd->qle;
            tle = swfwd->tle;
            gtle = swfwd->gtle;
            gscore = swfwd->gscore;
            
            // check whether we prefer to reach the end of the query
            if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
                a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
                a->truesc = a->score;
            } else { // to-end extension
                a->qb = 0, a->rb = s->rbeg - gtle;
                a->truesc = gscore;
            }
        }
        else{
            a->score = a->truesc = s->len * opt->a;
            a->qb = 0;
            a->rb = s->rbeg;
        }
        if(bwd->qlen!=0)//right extension
        {
            int qle, tle, qe, re, gtle, gscore, sc0 = swfwd->score;
            qe = s->qbeg + s->len;
            re = s->rbeg + s->len - rmax[0];
            a->score = swbwd->score;
            qle = swbwd->qle;
            tle = swbwd->tle;
            gtle = swbwd->gtle;
            gscore = swbwd->gscore;
            // similar to the above
            if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
                a->qe = qe + qle, a->re = rmax[0] + re + tle;
                a->truesc += a->score - sc0;
            } else { // to-end extension
                a->qe = l_query, a->re = rmax[0] + re + gtle;
                a->truesc += gscore - sc0;
            }
        }
        else{
            a->qe = l_query, a->re = s->rbeg + s->len;
        }
        // compute seedcov
        for (i = 0, a->seedcov = 0; i < c->n; ++i) {
            const mem_seed_t *t = &c->seeds[i];
            if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
                a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
        }
     //   a->w = aw;
        a->seedlen0 = s->len;
        a->frac_rep = c->frac_rep;
        if (bwa_verbose >= 4) printf("*** Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}\n", a->qb, a->qe, (long)a->rb, (long)a->re, a->score);
    }
}
#define min(a,b) a<b?a:b
#define max(a,b) a>b?a:b
void mem_chain2aln_filter(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_t *av_firstpass, mem_alnreg_v *av)
{
    int i, k;
    const mem_seed_t *s;
    //second pass
    // NEO:
    // external sorting
    // Generate str: str is an index,   high 32 bit is SW score (sorted)
    //                                  low 32 bit is real index
    uint64_t *srt;
    srt = malloc(c->n * 8);
    for (i = 0; i < c->n; ++i)
        srt[i] = (uint64_t)c->seeds[i].score<<32 | i;
    ks_introsort_64(c->n, srt);// NEO: srt in decending order //first by SW then by index
    
    for (k = c->n - 1; k >= 0; --k) {
        mem_alnreg_t *a;
        s = &c->seeds[(uint32_t)srt[k]];
        mem_alnreg_t *a_pre = &av_firstpass[(uint32_t)srt[k]];
        //mem_alnreg_t *a_prev = av_fistpass->a[(uint32_t)srt[k]];
        
        // NEO: this part is belong to CPU, should migrate this to the end of this function.
        // NEO:
        // should know how many seed would be drop in this place
        // Test if the seed is in future align
        for (i = 0; i < av->n; ++i) { // test whether extension has been made before
            mem_alnreg_t *p = &av->a[i];
            int64_t rd;
            int qd, w, max_gap, max_overlap;
            if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue; // not fully contained
            if (s->len - p->seedlen0 > .1 * l_query) continue; // this seed may give a better alignment
            // qd: distance ahead of the seed on query; rd: on reference
            qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
            max_overlap = min(qd,rd);
            assert(qd>=0);
            max_gap = cal_max_gap(opt, max_overlap); // the maximal gap allowed in regions ahead of the seed
          //  fprintf(stderr, "1 qd %d, rd %d, max_gap %d\n",qd,rd,max_gap);
            w = max_gap -max_overlap;//< p->w? max_gap : p->w; // bounded by the band width
            if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
            // similar to the previous four lines, but this time we look at the region behind
            qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
            max_overlap = min(qd,rd);
            max_gap = cal_max_gap(opt, max_overlap);
          //  fprintf(stderr, "2 qd %d, rd %d, max_gap %d\n",qd,rd,max_gap);
            w = max_gap-max_overlap;// < p->w? max_gap : p->w;
            if (qd - rd < w && rd - qd < w) break;
        }
   //     fprintf(stderr, "break is %d/%d", i<av->n,i);
        // NEO:
        // rescue the seed marked as overlap, if it would lead to a different result
        if (i < av->n) { // the seed is (almost) contained in an existing alignment; further testing is needed to confirm it is not leading to a different aln
            if (bwa_verbose >= 4)
                printf("** Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment [%d,%d) <=> [%ld,%ld)\n",
                       k, (long)s->len, (long)s->qbeg, (long)s->rbeg, av->a[i].qb, av->a[i].qe, (long)av->a[i].rb, (long)av->a[i].re);
            
            //NEO: block structure
            for (i = k + 1; i < c->n; ++i) { // check overlapping seeds in the same chain
                const mem_seed_t *t;
                if (srt[i] == 0) continue;
                t = &c->seeds[(uint32_t)srt[i]];
                if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODO: more efficient by early stopping
                if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
                if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
            }
            
            
            if (i == c->n) { // no overlapping seeds; then skip extension
                srt[k] = 0; // mark that seed extension has not been performed
                continue;
            }
            if (bwa_verbose >= 4)
                printf("** Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", k);
        }
        
        a = kv_pushp(mem_alnreg_t, *av);
        memcpy(a, a_pre, sizeof(mem_alnreg_t));
        
    }
    free(srt);
}

//this function should be change to batch and be seperatedinto three part, whith the mem_chain2aln needed to be changed into SIMD
mem_chain_v mem_gen_chains(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, void *buf)

{
    int i;
    mem_chain_v chn;
    for (i = 0; i < l_seq; ++i) // convert to 2-bit encoding if we have not done so
        seq[i] = seq[i] < 4? seq[i] : nst_nt4_table[(int)seq[i]];
    
    chn = mem_chain(opt, bwt, bns, l_seq, (uint8_t*)seq, buf);
    chn.n = mem_chain_flt(opt, chn.n, chn.a);
    mem_flt_chained_seeds(opt, bns, pac, l_seq, (uint8_t*)seq, chn.n, chn.a);
    if (bwa_verbose >= 4) mem_print_chain(bns, &chn);
    return chn;
}
//sw related to a query
typedef struct
{
    uint32_t l_query;//32
    size_t chn_count;//64
    //@para chn_count: a query's chain count;
    //@para max_id: a query's valid seed count;
    
    int64_t* g_maxs;
    
    uint8_t** rseqs;
    uint8_t** rseq_revs;
    
    const uint8_t *query;
    uint8_t* query_rev;
    
    mem_alnreg_t* g_av_firstpass;
    
    swseq_t* g_forward;
    swseq_t* g_backward;
    
    swrst_t* g_swfwd;
    swrst_t* g_swbwd;
    
    size_t* index;//save ist chain related value in index[i] location of g_forward/g_backward/g_swfwd/g_swbwd
    swrst_t* nxt_ext;
}qext_t;

void mem_chains2aln_init(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, mem_chain_v chn, qext_t* ext_val)
{
    int i;
    mem_alnreg_v regs;
    kv_init(regs);
    
    ext_val->l_query = l_seq;
    ext_val->query = (uint8_t*)seq;
    
    const uint8_t *query = (uint8_t*)seq;
    int l_query = l_seq;
    
    ext_val->chn_count = chn.n;
    
    //rev query init
    uint8_t *query_rev= 0;
    query_rev= malloc(l_query);
    for (i = 0; i < l_query; ++i) query_rev[i] = query[l_query - 1 - i];
    
    ext_val->query_rev=query_rev;
    
    //rmaxs
    int64_t * rmaxs = malloc(sizeof(int64_t)*2*chn.n);
    for(i=0; i<chn.n; i++)
    {
        int64_t l_pac = bns->l_pac;
        mem_chain2maxspan(opt, &chn.a[i],  l_query, l_pac, &rmaxs[2*i]);
    }
    ext_val->g_maxs = rmaxs;
    
    //reference
    uint8_t **rseqs = malloc(sizeof(uint8_t*)*chn.n);
    uint8_t **rseq_revs = malloc(sizeof(uint8_t*)*chn.n);
    for(int i=0; i<chn.n; i++)
    {
        int rid;
        int64_t rmax[2];
        rmax[0]=rmaxs[2*i];
        rmax[1]=rmaxs[2*i+1];
        rseqs[i] = bns_fetch_seq(bns, pac, &rmax[0], chn.a[i].seeds[0].rbeg, &rmax[1], &rid);
        assert(chn.a[i].rid == rid);
        uint32_t rseq_len = (uint32_t)(rmax[1]-rmax[0]);
        rseq_revs[i] = malloc(rseq_len);
        for (int k = 0; k < rseq_len; ++k) rseq_revs[i][k] = rseqs[i][rseq_len - 1 - k];
    }
    ext_val->rseqs=rseqs;
    ext_val->rseq_revs = rseq_revs;
    //scan
    size_t *index = malloc(sizeof(size_t)*(chn.n+1));
    index[0]=0;
    for(int i=1; i<=chn.n; i++)
    {
        const mem_chain_t *c = &chn.a[i-1];
        size_t tmp_pre = index[i-1];
        size_t tmp_n = c->n;
        tmp_pre+=tmp_n;
        index[i]=tmp_pre;
    }
    
    ext_val->index=index;
    
    size_t max_id =index[chn.n];
    swseq_t * g_forward = malloc(sizeof(swseq_t)*max_id);
    swseq_t * g_backward = malloc(sizeof(swseq_t)*max_id);
    swrst_t * g_swfwd = malloc(sizeof(swrst_t)*max_id);
    swrst_t * g_swbwd = malloc(sizeof(swrst_t)*max_id);
    mem_alnreg_t *g_av_firstpass = malloc(sizeof(mem_alnreg_t)*max_id);

    ext_val->g_swfwd = g_swfwd;
    ext_val->nxt_ext = g_swfwd;
    ext_val->g_swbwd = g_swbwd;
    ext_val->g_forward=g_forward;
    ext_val->g_backward=g_backward;
    ext_val->g_av_firstpass=g_av_firstpass;
    
    //main loop
    //SW init
    for (i = 0; i < chn.n; ++i) {
        const mem_chain_t *c = &chn.a[i];
        if (c->n == 0) continue;
        
        int64_t rmax[2];
        uint8_t *rseq = 0;
        uint8_t *rseq_rev = 0;
        
        swseq_t * forward = &g_forward[index[i]];
        swseq_t * backward = &g_backward[index[i]];
        swrst_t * swfwd = &g_swfwd[index[i]];
        swrst_t * swbwd = &g_swbwd[index[i]];
        
        rmax[0]=rmaxs[2*i];
        rmax[1]=rmaxs[2*i+1];
        
        rseq = rseqs[i];
        rseq_rev = rseq_revs[i];
        mem_chain2aln_preextent(opt, bns, pac, l_query, query, query_rev, rmax, rseq, rseq_rev, c, forward, backward, swfwd, swbwd);
    }
}

void mem_chains2aln_finalize(mem_chain_v chn,qext_t* ext_val)
{
    size_t *index = ext_val->index;

    int64_t * rmaxs = ext_val->g_maxs;
    mem_alnreg_t *g_av_firstpass = ext_val->g_av_firstpass;
    swseq_t * g_forward = ext_val->g_forward;
    swseq_t * g_backward = ext_val->g_backward;
    uint8_t **rseqs = ext_val->rseqs;
    uint8_t **rseq_revs = ext_val->rseq_revs;
    uint8_t *query_rev=ext_val->query_rev;
    for (int i = 0; i < chn.n; ++i) {
        free(chn.a[i].seeds);
    }
    free(chn.a);
    free(g_forward);
    free(g_backward);
    free(ext_val->g_swfwd);
    free(ext_val->g_swbwd);
    free(g_av_firstpass);
    free(index);
    for(int i=0; i<chn.n; i++)
    {
        free(rseqs[i]);
        free(rseq_revs[i]);
    }
    free(rseqs);
    free(rseq_revs);
    
    free(query_rev);
    free(rmaxs);
    
}
void mem_chains2aln_postextent(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, mem_chain_v chn,qext_t* ext_val ,mem_alnreg_v *av )
{
    int l_query=ext_val->l_query;
    const uint8_t *query = ext_val->query;
    mem_alnreg_t *g_av_firstpass = ext_val->g_av_firstpass;
    size_t *index = ext_val->index;
    int64_t * rmaxs = ext_val->g_maxs;
    //post SW
    for (int i=0; i<chn.n; ++i) {
        const mem_chain_t *c = &chn.a[i];
        if (c->n == 0) continue;
        mem_alnreg_t *av_firstpass =&g_av_firstpass[index[i]];
        swrst_t * swfwd = &ext_val->g_swfwd[index[i]];
        swrst_t * swbwd = &ext_val->g_swbwd[index[i]];
        mem_chain2aln_post(opt, bns, pac, l_query, query, c, &rmaxs[2*i], av_firstpass, swfwd, swbwd);
        mem_chain2aln_filter(opt, bns, pac, l_query, query, c, av_firstpass, av);
    }
}

void mem_chains2aln_sw(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, mem_chain_v chn, qext_t* ext_val)
{
    size_t *index = ext_val->index;
    int l_query=ext_val->l_query;
    const uint8_t *query = ext_val->query;
    
    //left extend
    for (int i = 0; i < chn.n; ++i) {
        const mem_chain_t *c = &chn.a[i];
        if (c->n == 0) continue;
        swrst_t * swfwd = &ext_val->g_swfwd[index[i]];
        
        //left extend
        mem_chain_extent(opt, bns, pac, l_query,query, c, swfwd);
    }
    for (int i = 0; i < chn.n; ++i) {
        const mem_chain_t *c = &chn.a[i];
        if (c->n == 0) continue;
        swrst_t * swfwd = &ext_val->g_swfwd[index[i]];
        swrst_t * swbwd = &ext_val->g_swbwd[index[i]];
        //init left extend
        for (int k = 0; k < c->n; ++k) {
            int pre_score= swfwd[k].score;
            swrst_t *swbwd_ = &swbwd[k];
            swbwd_->h0=pre_score;
        }
    }
    //right extent
    for (int i = 0; i < chn.n; ++i) {
        const mem_chain_t *c = &chn.a[i];
        if (c->n == 0) continue;
        swrst_t * swbwd = &ext_val->g_swbwd[index[i]];
        //right extend
        mem_chain_extent(opt, bns, pac, l_query,query, c, swbwd);
    }
}

mem_alnreg_v mem_chains2aln(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, mem_chain_v chn)
{
    //initialize
    mem_alnreg_v regs;
    kv_init(regs);
    
    qext_t* ext_val = malloc(sizeof(qext_t));
    
    mem_chains2aln_init(opt, bwt, bns, pac, l_seq, seq,chn, ext_val);
    
    //SW computation [GPU parallel]
    mem_chains2aln_sw(opt, bwt, bns, pac, l_seq, seq, chn, ext_val);
    //post SW
    mem_chains2aln_postextent(opt, bwt, bns, pac, l_seq, seq, chn, ext_val , &regs );
    
    //finalize
    mem_chains2aln_finalize( chn, ext_val);
    free(ext_val);
    return regs;
}



mem_alnreg_v mem_aln2regs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq, mem_alnreg_v regs)
{
    int i;
    regs.n = mem_sort_dedup_patch(opt, bns, pac, (uint8_t*)seq, regs.n, regs.a);
    if (bwa_verbose >= 4) {
        err_printf("* %ld chains remain after removing duplicated chains\n", regs.n);
        for (i = 0; i < regs.n; ++i) {
            mem_alnreg_t *p = &regs.a[i];
            printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
        }
    }
    for (i = 0; i < regs.n; ++i) {
        mem_alnreg_t *p = &regs.a[i];
        if (p->rid >= 0 && bns->anns[p->rid].is_alt)
            p->is_alt = 1;
    }
    return regs;
}


typedef struct {
    const mem_opt_t *opt;
    const bwt_t *bwt;
    const bntseq_t *bns;
    const uint8_t *pac;
    const mem_pestat_t *pes;
    smem_aux_t **aux;
    bseq1_t *seqs;
    mem_alnreg_v *regs;
  //  mem_chain_v *chn;//to do set to thread private
    
   // qext_t* ext_val;// to do set to thread private
    
    int64_t n_processed;
} worker_t_mod;
//NEO: this function need to be change to batch awared mode.


static void worker_aln2regs(void *data, int i, int tid)
{
    worker_t_mod *w = (worker_t_mod*)data;
     w->regs[i] = mem_aln2regs(w->opt, w->bwt, w->bns, w->pac, w->seqs[i].l_seq, w->seqs[i].seq, w->regs[i]);
}



static swrst_t* fwd(qext_t* ext_val,size_t id)
{
    return &ext_val->g_swfwd[id];
}

static swrst_t* bwd(qext_t* ext_val,size_t id)
{
    return &ext_val->g_swbwd[id];
}
//convert to SIMD
//init ref in db and query to qprof
//main process
//process result


//future SIMD
void mem_chain_extent_batch3(const mem_opt_t *opt, qext_t* ext_base, size_t* chn_idx, int batch, swrst_t*(*getItem)(qext_t*,size_t), int start,  int tid)
{
    //init
    size_t* batch_id = malloc((batch+1)*sizeof(size_t));
    batch_id[0]=0;
    size_t pre_id=0;
    for(int i=1; i<=batch; i++)
    {
        
        qext_t* ext_val = ext_base+i-1;
        size_t *ext_index = ext_val->index;
        size_t ext_size = ext_index[chn_idx[i-1]];// a seed related swrst count
        size_t cur_id = ext_size + pre_id;
        
        batch_id[i] = cur_id;
        pre_id = cur_id;
    }
    swrst_t* g_srt = malloc(batch_id[batch]*sizeof(swrst_t));
    for(int i=0; i<batch; i++)
    {
        qext_t* ext_val = ext_base+i;
        size_t *ext_index = ext_val->index;
        size_t ext_size = ext_index[chn_idx[i]];
        
        size_t ptr = batch_id[i];
        size_t idff = batch_id[i+1]-batch_id[i];
        assert(idff == ext_size);
        
        swrst_t * sws =getItem(ext_val,0);
        memcpy(g_srt+ptr, sws, ext_size*sizeof(swrst_t));
    }
    //SW
#ifdef DEBUG_SW
    kstring_t str={0,0,0};
    
    if(tid == 0&&start == 8000)
    {
        kputs("sw_start_",&str);
        kputl(start,&str);
        kputc('_',&str);
        kputl(tid,&str);
        kputc('_',&str);
        kputl(batch,&str);
        kputs(".bin",&str);
        store(g_srt,batch_id[batch],str.s);
    }
#endif
    
    ksw_extend_batch2(g_srt, (uint32_t)batch_id[batch], 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, opt->zdrop);
    
#ifdef DEBUG_SW
    if(tid == 0&&start == 8000)
    {
        free(str.s);
        str.l=0;
        str.m=0;
        str.s=0;
        kputs("sw_end_",&str);
        kputl(start,&str);
        kputc('_',&str);
        kputl(tid,&str);
        kputc('_',&str);
        kputl(batch,&str);
        kputs(".bin",&str);
        store(g_srt,batch_id[batch],str.s);
        free(str.s);
    }
#endif
    //finalize
    for(int i=0; i<batch; i++)
    {
        qext_t* ext_val = ext_base+i;
        size_t *ext_index = ext_val->index;
        size_t ext_size = ext_index[chn_idx[i]];
        
        size_t ptr = batch_id[i];
        
        swrst_t * sws =getItem(ext_val,0);
        memcpy(sws, g_srt+ptr, ext_size*sizeof(swrst_t));
    }
    
    free(batch_id);
    free(g_srt);
}

void mem_chain_extent_batch(const mem_opt_t *opt, qext_t* ext_base, size_t* chn_idx, int batch, swrst_t*(*getItem)(qext_t*,size_t), int start,  int tid)
{
    //init
    size_t* batch_id = malloc((batch+1)*sizeof(size_t));
    batch_id[0]=0;
    size_t pre_id=0;
    for(int i=1; i<=batch; i++)
    {
        qext_t* ext_val = ext_base+i-1;
        size_t *ext_index = ext_val->index;
        size_t ext_size = ext_index[chn_idx[i-1]];// a seed related swrst count
        size_t cur_id = ext_size + pre_id;
        batch_id[i] = cur_id;
        pre_id = cur_id;
    }
    swrst_t* g_srt = malloc(batch_id[batch]*sizeof(swrst_t));
    for(int i=0; i<batch; i++)
    {
        qext_t* ext_val = ext_base+i;
        size_t *ext_index = ext_val->index;
        size_t ext_size = ext_index[chn_idx[i]];
        
        size_t ptr = batch_id[i];
        size_t idff = batch_id[i+1]-batch_id[i];
        assert(idff == ext_size);
        
        swrst_t * sws =getItem(ext_val,0);
        memcpy(g_srt+ptr, sws, ext_size*sizeof(swrst_t));
    }
    //SW

    ksw_extend_batch2(g_srt, (uint32_t)batch_id[batch], 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, opt->zdrop);

    //finalize
    for(int i=0; i<batch; i++)
    {
        qext_t* ext_val = ext_base+i;
        size_t *ext_index = ext_val->index;
        size_t ext_size = ext_index[chn_idx[i]];
        
        size_t ptr = batch_id[i];
        
        swrst_t * sws =getItem(ext_val,0);
        memcpy(sws, g_srt+ptr, ext_size*sizeof(swrst_t));
    }
    
    free(batch_id);
    free(g_srt);
}

static void worker_mod_batch2(void *data, int start, int batch, int tid)
{
    worker_t_mod *w = (worker_t_mod*)data;
    qext_t* ext_base = malloc(sizeof(qext_t)*batch);//&w->ext_val[batch*tid];
    //gen chains
    size_t * chn_idx = malloc(sizeof(size_t)*batch);
    mem_chain_v * local_chn = malloc(sizeof(mem_chain_v)*batch);
    for(int i=start,j=0; j<batch; ++j,++i)
    {
        mem_chain_v chn =mem_gen_chains(w->opt, w->bwt, w->bns, w->pac, w->seqs[i].l_seq, w->seqs[i].seq, w->aux[tid]);
        chn_idx[j] = chn.n;
        local_chn[j]= chn;
    }
    //init sw
    for(int i=start,j=0; j<batch; ++j,++i)
    {
        qext_t* ext_val = ext_base+j;// &w->ext_val[batch*tid+j];
        
        mem_chains2aln_init(w->opt, w->bwt, w->bns, w->pac, w->seqs[i].l_seq, w->seqs[i].seq, local_chn[j], ext_val);
        
    }
    const mem_opt_t *opt = w->opt;
    
    //left extension
    mem_chain_extent_batch(opt, ext_base, chn_idx, batch, fwd,start, tid);
    //init right
    for(int j=0; j<batch; ++j)
    {
        qext_t* ext_val = ext_base+j;// &w->ext_val[batch*tid+j];
        
        size_t *index = ext_val->index;
        size_t max_idx = index[chn_idx[j]];
        swrst_t * swfwd = &ext_val->g_swfwd[0];
        swrst_t * swbwd = &ext_val->g_swbwd[0];
        for(int i=0; i<max_idx; i++)
        {
            int pre_score= swfwd[i].score;
            swrst_t *swbwd_ = &swbwd[i];
            swbwd_->h0=pre_score;
        }
    }
    //right extension
    mem_chain_extent_batch(opt, ext_base, chn_idx, batch, bwd,start,tid);
    //finalize
    for(int i=start,j=0; j<batch; ++j,++i)
    {
        qext_t* ext_val = ext_base+j;// &w->ext_val[batch*tid+j];
        
        mem_alnreg_v regs;
        kv_init(regs);
        mem_chains2aln_postextent(w->opt, w->bwt, w->bns, w->pac, w->seqs[i].l_seq, w->seqs[i].seq, local_chn[j], ext_val , &regs );
        //finalize
        w->regs[i] = regs;
    }
    for(int i=start,j=0; j<batch; ++j,++i)
    {
        qext_t* ext_val = ext_base+j;
        mem_chains2aln_finalize( local_chn[j], ext_val);
    }
    free(ext_base);
    free(local_chn);
    free(chn_idx);
}
typedef struct
{
    mem_alnreg_v*av;
    mem_alnreg_t*a;
    mem_seed_t* s;
    
    const mem_chain_t*c;
    int64_t *rmax;
    
    uint8_t *rseq;
    const uint8_t *query ;
    
    int l_query;
    
}sw_itv_val;
typedef struct
{
    size_t n, m;
    sw_itv_val* a;
}sw_itv_vec;

static void worker_mod_batch(void *data, int start, int batch, int tid)
{
    worker_t_mod *w = (worker_t_mod*)data;
    mem_chain_v *local_chn = malloc(sizeof(mem_chain_v)*batch);
    mem_alnreg_v *local_regs = malloc(sizeof(mem_alnreg_v)*batch);
    
    const mem_opt_t *opt = w->opt;
    const bntseq_t *bns = w->bns;
    const uint8_t *pac = w->pac;
    size_t * global_chn_id = malloc(sizeof(size_t)*(batch+1));
    global_chn_id[0]=0;
    for(int i=start, j=0; j<batch; j++,i++)
    {
        mem_chain_v tmp_chn = mem_gen_chains(w->opt, w->bwt, w->bns, w->pac, w->seqs[i].l_seq, w->seqs[i].seq, w->aux[tid]);
        local_chn[j] = tmp_chn;
        global_chn_id[j+1] = global_chn_id[j]+tmp_chn.n;
    }
    
    //initialize of SW
    for(int i=start, j=0; j<batch; j++,i++)
    {
        mem_alnreg_v* regs = &local_regs[j];
        kv_init(*regs);
    }
    //batch values:

    int64_t* g_rmaxs = malloc(sizeof(int64_t)*2*global_chn_id[batch]);
    mem_alnreg_v* global_regs = malloc(sizeof(mem_alnreg_v)*global_chn_id[batch]);
    mem_chain_t* global_chain_t = malloc(sizeof(mem_chain_t)*global_chn_id[batch]);
//    swseq_t* g_forward = malloc(sizeof(swseq_t*)*global_chn_id[batch]);
//    swseq_t* g_backward = malloc(sizeof(swseq_t*)*global_chn_id[batch]);
//    
//    swrst_t* g_swfwd = malloc(sizeof(swrst_t*)*global_chn_id[batch]);
//    swrst_t* g_swbwd = malloc(sizeof(swrst_t*)*global_chn_id[batch]);
//    
 //   swseq_t * batch_swseq = malloc(sizeof(swseq_t)*batch);
//    swrst_t * batch_swfwd = malloc(sizeof(swrst_t)*batch);
    for(int i=0; i<global_chn_id[batch]; i++)
    {
        kv_init(global_regs[i]);
    }
 //   uint64_t **global_srt = malloc(sizeof(uint64_t*)*batch);
   // uint8_t **global_rseq = malloc(sizeof(uint8_t*)*batch);
    int* global_seqlen = malloc(sizeof(int)*global_chn_id[batch]);
    char* * global_seq = malloc(sizeof(char*)*global_chn_id[batch]);
    uint8_t **global_rseq = malloc(sizeof(uint8_t*)*global_chn_id[batch]);
    uint64_t ** global_srt = malloc(sizeof(uint64_t*)* global_chn_id[batch]);
    //expend the second loop
    for(int i=start, j=0; j<batch; j++,i++)
    {
        mem_chain_v chn_v = local_chn[j];
        for (int l_chn_id = 0; l_chn_id < chn_v.n; ++l_chn_id) {
            
            global_chain_t[(global_chn_id[j]+l_chn_id)] = chn_v.a[l_chn_id];
            global_seqlen[(global_chn_id[j]+l_chn_id)] =w->seqs[i].l_seq;
            global_seq[(global_chn_id[j]+l_chn_id)]=w->seqs[i].seq;
        }
    }
    
    //for(int j=0; j<batch; j++)
    int SW_batch = batch*2;
    int max_c = global_chn_id[batch];
    int seg = (max_c+SW_batch-1)/SW_batch;
    int next_process = SW_batch;
    int left = max_c;
    //for(int g_c_id=0; g_c_id<max_c; g_c_id++)
    
    int *seeds_idx = malloc(sizeof(int)*SW_batch);
    sw_itv_vec sw_nxt_process;
    kv_init(sw_nxt_process);
  //  int *seeds_end = malloc(sizeof(int)*SW_batch);
    
    for(int seg_id=0; seg_id<seg; seg_id++)
    {
        next_process = next_process<left?next_process:left;
        left-=SW_batch;
        for(int i=0; i<next_process; i++)//process batch of data
        {
            {
                int g_c_id = SW_batch*seg_id+i;
                mem_chain_t*p = &global_chain_t[g_c_id];//&chn_v.a[l_chn_id];
                
                int l_query = global_seqlen[g_c_id];//l_seq;
               
                const mem_chain_t*c = p;
                
                int i, rid;// aw[2]; // aw: actual bandwidth used in extension
                int64_t l_pac = bns->l_pac,/* rmax[2],*/  max = 0;
                
                
                int64_t *rmax = &g_rmaxs[(g_c_id)*2];//&b_rmaxs[j*2];
                
                if (c->n == 0) return;
                // get the max possible span
                // NEO: should set on CPU
                // NEO: can we orgnize query by max possible span?
                rmax[0] = l_pac<<1; rmax[1] = 0;
                for (i = 0; i < c->n; ++i) {
                    int64_t b, e;
                    const mem_seed_t *t = &c->seeds[i];
                    b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
                    e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
                    rmax[0] = rmax[0] < b? rmax[0] : b;
                    rmax[1] = rmax[1] > e? rmax[1] : e;
                    if (t->len > max) max = t->len;
                }
                rmax[0] = rmax[0] > 0? rmax[0] : 0;
                rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
                if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
                    if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
                    else rmax[0] = l_pac;
                }
                // retrieve the reference sequence
                uint8_t *rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);//NEO: potentially OOM, in every 10MB batch, average 67MB
                global_rseq[g_c_id] = rseq;
                // global_rseq[j]=rseq;
                assert(c->rid == rid);
                
                // NEO:
                // external sorting
                // Generate str: str is an index,   high 32 bit is SW score (sorted)
                //                                  low 32 bit is real index
                uint64_t * srt = malloc(c->n * 8);
                   //     global_srt[j] = srt;
                for (i = 0; i < c->n; ++i)
                    srt[i] = (uint64_t)c->seeds[i].score<<32 | i;
                
                ks_introsort_64(c->n, srt);// NEO: srt in decending order
                global_srt[g_c_id] = srt;
            }
        }
       
        /*
         NEO's plan:
         outer loop: while(1)
         inner loop, filter produce next seeds,
         sizeof(seeds)==0 break;
         extend base on seeds.
         post process.
         */
        for(int cur_process_id=0; cur_process_id<next_process; cur_process_id++)//process batch of data
        {
            int g_c_id = SW_batch*seg_id+cur_process_id;
            seeds_idx[cur_process_id] = global_chain_t[g_c_id].n-1;
        }
        while(1)
        {
            for(int cur_process_id=0; cur_process_id<next_process; cur_process_id++)//process batch of data
            {
            // NEO: should do modification in this part in the future

                int g_c_id = SW_batch*seg_id+cur_process_id;
                uint64_t *srt = global_srt[g_c_id];
                uint8_t *rseq = global_rseq[g_c_id];
                int64_t *rmax = &g_rmaxs[(g_c_id)*2];//&b_rmaxs[j*2];
                const uint8_t *query = (uint8_t*)global_seq[g_c_id];//seq;
               // mem_alnreg_v* regs = ;
                mem_alnreg_v*av = &global_regs[g_c_id];
                mem_chain_t*p = &global_chain_t[g_c_id];//&chn_v.a[l_chn_id];
                const mem_chain_t*c = p;
                int l_query = global_seqlen[g_c_id];//l_seq;
                
                
                for (; seeds_idx[cur_process_id] >= 0; --seeds_idx[cur_process_id]) {
                    {
                        mem_seed_t *s = &c->seeds[(uint32_t)srt[seeds_idx[cur_process_id]]];
                        // NEO: this part is belong to CPU, should migrate this to the end of this function.
                        // Test if the seed is in future align
                        // NEO: Filter
                        int i;
                        for (i = 0; i < av->n; ++i) { // test whether extension has been made before
                                    mem_alnreg_t *p = &av->a[i];
                                    int64_t rd;
                                    int qd, w, max_gap, max_overlap;
                                    if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue; // not fully contained
                                    if (s->len - p->seedlen0 > .1 * l_query) continue; // this seed may give a better alignment
                                    // qd: distance ahead of the seed on query; rd: on reference
                                    qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
                                    max_overlap = min(qd,rd);
                                    max_gap = cal_max_gap(opt, max_overlap); // the maximal gap allowed in regions ahead of the seed
                                    w = max_gap < p->w? max_gap : p->w;
                                    //w = max_gap -max_overlap;//< p->w? max_gap : p->w; // bounded by the band width
                                    if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
                                    // similar to the previous four lines, but this time we look at the region behind
                                    qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
                                    max_overlap = min(qd,rd);
                                    max_gap = cal_max_gap(opt, max_overlap);
                                    w = max_gap < p->w? max_gap : p->w;
                                    //w = max_gap-max_overlap;// < p->w? max_gap : p->w;
                                    if (qd - rd < w && rd - qd < w) break;
                                }
                                
                        // NEO:
                        // rescue the seed marked as overlap, if it would lead to a different result
                        if (i < av->n) { // the seed is (almost) contained in an existing alignment; further testing is needed to confirm it is not leading to a different aln
                                    if (bwa_verbose >= 4)
                                        printf("** Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment [%d,%d) <=> [%ld,%ld)\n",
                                               seeds_idx[cur_process_id], (long)s->len, (long)s->qbeg, (long)s->rbeg, av->a[i].qb, av->a[i].qe, (long)av->a[i].rb, (long)av->a[i].re);
                                    
                                    //NEO: block structure
                                    for (i = seeds_idx[cur_process_id] + 1; i < c->n; ++i) { // check overlapping seeds in the same chain
                                        const mem_seed_t *t;
                                        if (srt[i] == 0) continue;
                                        t = &c->seeds[(uint32_t)srt[i]];
                                        if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODO: more efficient by early stopping
                                        if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
                                        if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
                                    }
                                    
                                    if (i == c->n) { // no overlapping seeds; then skip extension
                                        srt[seeds_idx[cur_process_id]] = 0; // mark that seed extension has not been performed
                                        continue;
                                    }
                                    if (bwa_verbose >= 4)
                                        printf("** Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", seeds_idx[cur_process_id]);
                                }
                        sw_itv_val tmp_sw_itv;
                        tmp_sw_itv.s = s;
                        tmp_sw_itv.av = av;
                        
//                        mem_alnreg_t* a = kv_pushp(mem_alnreg_t, *av);
//                        memset(a, 0, sizeof(mem_alnreg_t));
//                        a->w  = opt->w;
//                        a->score = a->truesc = -1;
//                        a->rid = c->rid;
//                        tmp_sw_itv.a=a;
                        
                        tmp_sw_itv.query = query;
                        tmp_sw_itv.rseq = rseq;
                        tmp_sw_itv.c = c;
                        tmp_sw_itv.l_query = l_query;
                        tmp_sw_itv.rmax=rmax;
                        //init the values used in SW extent
                        kv_push(sw_itv_val, sw_nxt_process, tmp_sw_itv);
                    }
                }
            }
                    /**********************/
#ifdef DEBUG
            if(sw_nxt_process.n<SW_batch)
                fprintf(stderr,"counting %ld\n",sw_nxt_process.n);
#endif
            if(sw_nxt_process.n==0)goto endwhile;
                    // mem_alnreg_v*av = regs;
            
            swseq_t* b_sw_seq_left = malloc(sizeof(swseq_t)*sw_nxt_process.n);
            swrst_t* b_sw_vals_left = malloc(sizeof(swrst_t)*sw_nxt_process.n);
            swseq_t* b_sw_seq_right = malloc(sizeof(swseq_t)*sw_nxt_process.n);
            swrst_t* b_sw_vals_right = malloc(sizeof(swrst_t)*sw_nxt_process.n);
            memset(b_sw_seq_left, 0, sizeof(swseq_t));
            memset(b_sw_seq_right, 0, sizeof(swseq_t));
            for(int cur_ptr=0; cur_ptr<sw_nxt_process.n; cur_ptr++)//init sw related values
            {
                swseq_t* cur_seq = &b_sw_seq_right[cur_ptr];
                swrst_t* cur_srt = &b_sw_vals_right[cur_ptr];
                cur_srt->sw_seq=cur_seq;
                cur_seq = &b_sw_seq_left[cur_ptr];
                cur_srt = &b_sw_vals_left[cur_ptr];
                cur_srt->sw_seq=cur_seq;
            }
            for(int cur_ptr=0; cur_ptr<sw_nxt_process.n; cur_ptr++) // left extension init
            {
                int64_t tmp;
                sw_itv_val tmp_sw_itv = sw_nxt_process.a[cur_ptr];
                mem_seed_t* s = tmp_sw_itv.s;
                const uint8_t *query =tmp_sw_itv.query;
                uint8_t *rseq=tmp_sw_itv.rseq;
                int64_t* rmax = tmp_sw_itv.rmax;
                
                swrst_t* cur_srt_l = &b_sw_vals_left[cur_ptr];
                swseq_t* cur_seq_l = cur_srt_l->sw_seq;
                
                // MAIN SW
                if (s->qbeg) { // left extension init
                    uint8_t *rs, *qs;
                    
                    qs = malloc(s->qbeg);
                    for (int i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];//query
                    tmp = s->rbeg - rmax[0];
                    rs = malloc(tmp);
                    for (int i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];//rseq
                    cur_seq_l->qlen = s->qbeg;
                    cur_seq_l->query = qs;
                    cur_seq_l->rlen = tmp;
                    cur_seq_l->ref = rs;
                    cur_srt_l->h0 = s->len * opt->a;
                }
                else{
                    cur_seq_l->qlen = 0;//NEO: just set a flag
                    cur_seq_l->rlen = 0;
                }
            }
//            ksw_extend_batch2(b_sw_vals_left, (uint32_t)sw_nxt_process.n, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, opt->zdrop);
            ksw_extend_batchw(b_sw_vals_left, (uint32_t)sw_nxt_process.n, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, opt->w, opt->pen_clip5, opt->zdrop);

            for(int cur_ptr=0; cur_ptr<sw_nxt_process.n; cur_ptr++)// right extention init
            {
                sw_itv_val tmp_sw_itv = sw_nxt_process.a[cur_ptr];
                mem_seed_t* s = tmp_sw_itv.s;
                const uint8_t *query =tmp_sw_itv.query;
                uint8_t *rseq=tmp_sw_itv.rseq;
                int l_query = tmp_sw_itv.l_query;
                int64_t* rmax = tmp_sw_itv.rmax;
                
                swrst_t* cur_srt_l = &b_sw_vals_left[cur_ptr];
                swseq_t* cur_seq_l = cur_srt_l->sw_seq;
                
                if (cur_seq_l->qlen==0) {
                    cur_srt_l->score = s->len * opt->a;
                }
                
                swrst_t* cur_srt_r = &b_sw_vals_right[cur_ptr];
                swseq_t* cur_seq_r = cur_srt_r->sw_seq;
                if (s->qbeg + s->len != l_query) { // right extension init
//                    int qle, tle, qe, re, gtle, gscore, sc0 = a->score;
                    int qe,re;
                    qe = s->qbeg + s->len;
                    re = s->rbeg + s->len - rmax[0];
                    assert(re >= 0);
                    //NEO: warp or block
                    cur_seq_r->qlen = l_query - qe;
                    cur_seq_r->query = query + qe;
                    cur_seq_r->rlen = rmax[1] - rmax[0] - re;
                    cur_seq_r->ref = rseq + re;
                    cur_srt_r->h0 = b_sw_vals_left[cur_ptr].score;
                }
                else{
                    cur_seq_r->qlen = 0;
                    cur_seq_r->rlen = 0;
                }
            }

//            ksw_extend_batch2(b_sw_vals_right, (uint32_t)sw_nxt_process.n, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, opt->zdrop);
            ksw_extend_batchw(b_sw_vals_right, (uint32_t)sw_nxt_process.n, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, opt->w, opt->pen_clip3, opt->zdrop);
            for(int cur_ptr=0; cur_ptr<sw_nxt_process.n; cur_ptr++)// post process (left & right)
            {
                int max_off[2];
                sw_itv_val tmp_sw_itv = sw_nxt_process.a[cur_ptr];
                mem_alnreg_v*av = tmp_sw_itv.av;
                mem_seed_t* s = tmp_sw_itv.s;
                const mem_chain_t* c = tmp_sw_itv.c;
                int l_query = tmp_sw_itv.l_query;
                int64_t* rmax = tmp_sw_itv.rmax;
                swrst_t* cur_srt_r = &b_sw_vals_right[cur_ptr];
                swseq_t* cur_seq_r = cur_srt_r->sw_seq;
                
                mem_alnreg_t* a = kv_pushp(mem_alnreg_t, *av);
                memset(a, 0, sizeof(mem_alnreg_t));
                a->w  = opt->w;
                a->score = a->truesc = -1;
                a->rid = c->rid;
                swrst_t* cur_srt_l = &b_sw_vals_left[cur_ptr];
                swseq_t* cur_seq_l = cur_srt_l->sw_seq;
                if (cur_seq_l->qlen!=0) {// left extension poster
                    a->score = cur_srt_l->score;
                    int qle, tle, gtle, gscore;
                    qle = cur_srt_l->qle;
                    tle = cur_srt_l->tle;
                    gtle = cur_srt_l->gtle;
                    gscore = cur_srt_l->gscore;
                    max_off[0]=cur_srt_l->max_off;
                    
                    // check whether we prefer to reach the end of the query
                    if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
                        a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
                        a->truesc = a->score;
                    } else { // to-end extension
                        a->qb = 0, a->rb = s->rbeg - gtle;
                        a->truesc = gscore;
                    }
                    uint8_t *rs, *qs;
                    rs = (uint8_t *)cur_seq_l->ref;
                    qs = (uint8_t *)cur_seq_l->query;
                    free(qs); free(rs);
                }
                else
                {
                    a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;
                    // swrst_t* cur_srt = &b_sw_vals_left[cur_ptr];
                }
                
                if(cur_seq_r->qlen!=0)//right extension post process
                {
                    int qle, tle, qe, re, gtle, gscore;
                    int sc0=cur_srt_r->h0;
                    a->score = cur_srt_r->score;
                    qle = cur_srt_r->qle;
                    tle = cur_srt_r->tle;
                    gtle = cur_srt_r->gtle;
                    gscore = cur_srt_r->gscore;
                    max_off[1]=cur_srt_r->max_off;
                    qe = s->qbeg + s->len;
                    re = s->rbeg + s->len - rmax[0];
                            // similar to the above
                    if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
                        a->qe = qe + qle, a->re = rmax[0] + re + tle;
                        a->truesc += a->score - sc0;
                    } else { // to-end extension
                        a->qe = l_query, a->re = rmax[0] + re + gtle;
                        a->truesc += gscore - sc0;
                    }
                }
                else
                {
                    a->qe = l_query, a->re = s->rbeg + s->len;
                }
                        // compute seedcov
                int i;
                for (i = 0, a->seedcov = 0; i < c->n; ++i) {
                    const mem_seed_t *t = &c->seeds[i];
                    if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
                        a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
                }
                a->seedlen0 = s->len;
                a->frac_rep = c->frac_rep;//c
            }
            
            sw_nxt_process.n=0;//set zero
            free(b_sw_seq_left);
            free(b_sw_vals_left);
            free(b_sw_seq_right);
            free(b_sw_vals_right);
        }
    endwhile:
        for(int i=0; i<next_process; i++)//process batch of data
        {
            {
                int g_c_id = SW_batch*seg_id+i;
                uint64_t *srt = global_srt[g_c_id];
                uint8_t *rseq = global_rseq[g_c_id];
                free(srt);
                free(rseq);
            }
        }
        
    }
    
    //save result
    //expanded regs (global_regs) to original regs (local_regs)
    for(int i=start, j=0; j<batch; j++,i++)
    {
        mem_alnreg_v* cur_regs = &local_regs[j];
        for (int g_chn_id = global_chn_id[j];  g_chn_id< global_chn_id[j+1]; ++g_chn_id) {
            mem_alnreg_v* regs = &global_regs[g_chn_id];
            for(int i=0; i<regs->n; i++)
            {
                kv_push(mem_alnreg_t, *cur_regs, regs->a[i]);
            }
            free(regs->a);
        }
    }
    free(global_regs);
  //  free(global_srt);


    for(int i=start, j=0; j<batch; j++,i++)
    {
        char *seq = w->seqs[i].seq;
        mem_alnreg_v* regs = local_regs+j;
        regs->n = mem_sort_dedup_patch(opt, bns, pac, (uint8_t*)seq, regs->n, regs->a);
        if (bwa_verbose >= 4) {
            err_printf("* %ld chains remain after removing duplicated chains\n", regs->n);
            for (int i = 0; i < regs->n; ++i) {
                mem_alnreg_t *p = &regs->a[i];
                printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
            }
        }
        for (int i = 0; i < regs->n; ++i) {
            mem_alnreg_t *p = &regs->a[i];
            if (p->rid >= 0 && bns->anns[p->rid].is_alt)
                p->is_alt = 1;
        }
        w->regs[i] = *regs;
    }
    
    //finalize
    for(int j=0; j<batch; j++)
    {
        mem_chain_v chn = local_chn[j];
        for (int i = 0; i < chn.n; ++i) {
            free(chn.a[i].seeds);
        }
        free(chn.a);
    }
    free(local_chn);
    free(local_regs);
 //   free(batch_swfwd);
 //   free(batch_swseq);
    free(global_chn_id);
    free(global_chain_t);
    free(global_rseq);
    free(global_seqlen);
    free(global_seq);
    free(g_rmaxs);
    free(global_srt);
    free(seeds_idx);
    free(sw_nxt_process.a);
   // free(seeds_end);
}

void chainging_batch(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *seqs, smem_aux_t*aux, int batch, mem_chain_v *local_chnvs)
{
    for(int i=0; i<batch; i++)
    {
        int l_seq = seqs[i].l_seq;
        char *seq = seqs[i].seq;
        mem_chain_v chnv;
        chnv= mem_gen_chains(opt, bwt, bns, pac, l_seq, seq, aux);
        local_chnvs[i] =chnv;
    }
}


void seed_extension_batch(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *seqs, smem_aux_t*aux, int batch, mem_chain_v *local_chnvs, mem_alnreg_v *local_regvs)
{

    for(int i=0; i<batch; i++)
    {
        mem_alnreg_v* p_regs = local_regvs+i;
        kv_init(*p_regs);
    }
    int64_t** read_rmaxs = malloc(sizeof(int64_t)*batch);
    
    //init rmaxs (max possible span)
    /*
     NEO:
     @para rmax[2] {thread private}:
     [0]: begin
     [1]: end
     */
    for(int batch_id=0; batch_id<batch; batch_id++)//read
    {
        int l_seq = seqs[batch_id].l_seq;
        mem_chain_v chnv = local_chnvs[batch_id];
        read_rmaxs[batch_id] = malloc(sizeof(int64_t)*chnv.n*2);
        int64_t* chnv_rmaxs = read_rmaxs[batch_id];
        for (int chain_id = 0; chain_id < chnv.n; ++chain_id) {//chain inside read
            mem_chain_t *p = &chnv.a[chain_id];
            int l_query = l_seq;
            const mem_chain_t*c = p;
            
            int64_t *rmax = chnv_rmaxs+2*chain_id;
            int i; // aw: actual bandwidth used in extension
            int64_t l_pac = bns->l_pac, max = 0;
            rmax[0] = l_pac<<1; rmax[1] = 0;
            for (i = 0; i < c->n; ++i) {
                int64_t b, e;
                const mem_seed_t *t = &c->seeds[i];
                b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
                e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
                rmax[0] = rmax[0] < b? rmax[0] : b;
                rmax[1] = rmax[1] > e? rmax[1] : e;
                if (t->len > max) max = t->len;
            }
            rmax[0] = rmax[0] > 0? rmax[0] : 0;
            rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
            if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
                if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
                else rmax[0] = l_pac;
            }
        }
        
    }
    
    
    for(int batch_id=0; batch_id<batch; batch_id++)//read
    {
        int l_seq = seqs[batch_id].l_seq;
        char *seq = seqs[batch_id].seq;
        
        mem_chain_v chnv = local_chnvs[batch_id];
        mem_alnreg_v* p_regs = local_regvs+batch_id;
        
        int64_t* chnv_rmaxs = read_rmaxs[batch_id];
        
        for (int chain_id = 0; chain_id < chnv.n; ++chain_id) {//chain inside read
            mem_chain_t *p = &chnv.a[chain_id];
            if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", chain_id);
            {
                const uint8_t * query = (uint8_t*)seq;
                int l_query = l_seq;
                const mem_chain_t*c = p;
                mem_alnreg_v*av = p_regs;
                
                int64_t *rmax = chnv_rmaxs+2*chain_id;
                int i, k, rid, max_off[2], aw[2]; // aw: actual bandwidth used in extension
   //             int64_t l_pac = bns->l_pac,  tmp, max = 0;
                int64_t tmp;
                const mem_seed_t *s;
                uint8_t *rseq = 0;
                uint64_t *srt;
                
                if (c->n == 0) return;
                // retrieve the reference sequence
                rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);//NEO: potentially OOM, in every 10MB batch, average 67MB
                assert(c->rid == rid);
                
                // NEO:
                // external sorting
                // Generate str: str is an index,   high 32 bit is SW score (sorted)
                //                                  low 32 bit is real index
                srt = malloc(c->n * 8);
                for (i = 0; i < c->n; ++i)
                    srt[i] = (uint64_t)c->seeds[i].score<<32 | i;
                ks_introsort_64(c->n, srt);// NEO: srt in decending order
                
                
                // NEO: should do modification in this part in the future
                for (k = c->n - 1; k >= 0; --k) {//seed inside chain
                    mem_alnreg_t *a;
                    s = &c->seeds[(uint32_t)srt[k]];
                    
                    // NEO: this part is belong to CPU, should migrate this to the end of this function.
                    // NEO:
                    // should know how many seed would be drop in this place
                    // Test if the seed is in future align
                    for (i = 0; i < av->n; ++i) { // test whether extension has been made before
                        mem_alnreg_t *p = &av->a[i];
                        int64_t rd;
                        int qd, w, max_gap;
                        if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue; // not fully contained
                        if (s->len - p->seedlen0 > .1 * l_query) continue; // this seed may give a better alignment
                        // qd: distance ahead of the seed on query; rd: on reference
                        qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
                        max_gap = cal_max_gap(opt, qd < rd? qd : rd); // the maximal gap allowed in regions ahead of the seed
                        w = max_gap < p->w? max_gap : p->w; // bounded by the band width
                        if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
                        // similar to the previous four lines, but this time we look at the region behind
                        qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
                        max_gap = cal_max_gap(opt, qd < rd? qd : rd);
                        w = max_gap < p->w? max_gap : p->w;
                        if (qd - rd < w && rd - qd < w) break;
                    }
                    // NEO:
                    // rescue the seed marked as overlap, if it would lead to a different result
                    if (i < av->n) { // the seed is (almost) contained in an existing alignment; further testing is needed to confirm it is not leading to a different aln
                        if (bwa_verbose >= 4)
                            printf("** Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment [%d,%d) <=> [%ld,%ld)\n",
                                   k, (long)s->len, (long)s->qbeg, (long)s->rbeg, av->a[i].qb, av->a[i].qe, (long)av->a[i].rb, (long)av->a[i].re);
                        
                        //NEO: block structure
                        for (i = k + 1; i < c->n; ++i) { // check overlapping seeds in the same chain
                            const mem_seed_t *t;
                            if (srt[i] == 0) continue;
                            t = &c->seeds[(uint32_t)srt[i]];
                            if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODO: more efficient by early stopping
                            if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
                            if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
                        }
                        
                        if (i == c->n) { // no overlapping seeds; then skip extension
                            srt[k] = 0; // mark that seed extension has not been performed
                            continue;
                        }
                        if (bwa_verbose >= 4)
                            printf("** Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", k);
                    }
                    
                    a = kv_pushp(mem_alnreg_t, *av);
                    memset(a, 0, sizeof(mem_alnreg_t));
                    a->w = aw[0] = aw[1] = opt->w;
                    a->score = a->truesc = -1;
                    a->rid = c->rid;
                    
                    if (bwa_verbose >= 4) err_printf("** ---> Extending from seed(%d) [%ld;%ld,%ld] @ %s <---\n", k, (long)s->len, (long)s->qbeg, (long)s->rbeg, bns->anns[c->rid].name);
                    if (s->qbeg) { // left extension
                        uint8_t *rs, *qs;
                        int qle, tle, gtle, gscore;
                        qs = malloc(s->qbeg);
                        for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
                        tmp = s->rbeg - rmax[0];
                        rs = malloc(tmp);
                        for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];
                        for (i = 0; i < MAX_BAND_TRY; ++i) {
                            int prev = a->score;
                            aw[0] = opt->w << i;
                            if (bwa_verbose >= 4) {
                                int j;
                                printf("*** Left ref:   "); for (j = 0; j < tmp; ++j) putchar("ACGTN"[(int)rs[j]]); putchar('\n');
                                printf("*** Left query: "); for (j = 0; j < s->qbeg; ++j) putchar("ACGTN"[(int)qs[j]]); putchar('\n');
                            }
                            //NEO: the most time consuming part
                            a->score = ksw_extend2(s->qbeg, qs, tmp, rs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, s->len * opt->a, &qle, &tle, &gtle, &gscore, &max_off[0]);
                            if (bwa_verbose >= 4) { printf("*** Left extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[0], max_off[0]); fflush(stdout); }
                            if (a->score == prev || max_off[0] < (aw[0]>>1) + (aw[0]>>2)) break;
                        }
                        // check whether we prefer to reach the end of the query
                        if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
                            a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
                            a->truesc = a->score;
                        } else { // to-end extension
                            a->qb = 0, a->rb = s->rbeg - gtle;
                            a->truesc = gscore;
                        }
                        free(qs); free(rs);
                    } else a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;
                    
                    if (s->qbeg + s->len != l_query) { // right extension
                        int qle, tle, qe, re, gtle, gscore, sc0 = a->score;
                        qe = s->qbeg + s->len;
                        re = s->rbeg + s->len - rmax[0];
                        assert(re >= 0);
                        
                        //NEO: warp or block
                        for (i = 0; i < MAX_BAND_TRY; ++i) {
                            int prev = a->score;
                            aw[1] = opt->w << i;
                            if (bwa_verbose >= 4) {
                                int j;
                                printf("*** Right ref:   "); for (j = 0; j < rmax[1] - rmax[0] - re; ++j) putchar("ACGTN"[(int)rseq[re+j]]); putchar('\n');
                                printf("*** Right query: "); for (j = 0; j < l_query - qe; ++j) putchar("ACGTN"[(int)query[qe+j]]); putchar('\n');
                            }
                            a->score = ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
                            if (bwa_verbose >= 4) { printf("*** Right extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[1], max_off[1]); fflush(stdout); }
                            if (a->score == prev || max_off[1] < (aw[1]>>1) + (aw[1]>>2)) break;
                        }
                        
                        // similar to the above
                        if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
                            a->qe = qe + qle, a->re = rmax[0] + re + tle;
                            a->truesc += a->score - sc0;
                        } else { // to-end extension
                            a->qe = l_query, a->re = rmax[0] + re + gtle;
                            a->truesc += gscore - sc0;
                        }
                    } else a->qe = l_query, a->re = s->rbeg + s->len;
                    if (bwa_verbose >= 4) printf("*** Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}_bandwidth={%d,%d}\n", a->qb, a->qe, (long)a->rb, (long)a->re, a->score, aw[0], aw[1]);
                    
                    // compute seedcov
                    for (i = 0, a->seedcov = 0; i < c->n; ++i) {
                        const mem_seed_t *t = &c->seeds[i];
                        if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
                            a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
                    }
                    a->w = aw[0] > aw[1]? aw[0] : aw[1];
                    a->seedlen0 = s->len;
                    
                    a->frac_rep = c->frac_rep;
                }
                free(srt); free(rseq);
            }

        }
      //  free(chnv_rmaxs);
    }
    
    //finalize rmaxs
    for(int i=0; i<batch; i++)//read
    {
        free(read_rmaxs[i]); 
    }
    
    free(read_rmaxs);
}

void post_extensiong_batch(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *seqs, int batch,  mem_alnreg_v *local_regvs,  mem_alnreg_v *global_regvs)
{
    for(int i=0; i<batch; i++)
    {
        char *seq = seqs[i].seq;
        mem_alnreg_v regs = local_regvs[i];
        regs.n = mem_sort_dedup_patch(opt, bns, pac, (uint8_t*)seq, regs.n, regs.a);
        if (bwa_verbose >= 4) {
            err_printf("* %ld chains remain after removing duplicated chains\n", regs.n);
            for (int i = 0; i < regs.n; ++i) {
                mem_alnreg_t *p = &regs.a[i];
                printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
            }
        }
        for (int i = 0; i < regs.n; ++i) {
            mem_alnreg_t *p = &regs.a[i];
            if (p->rid >= 0 && bns->anns[p->rid].is_alt)
                p->is_alt = 1;
        }
        global_regvs[i] = regs;
    }
}

static void worker1_batch(void *data, int start, int batch, int tid)
{
    worker_t_mod *w = (worker_t_mod*)data;
    
    mem_chain_v *local_chnvs = malloc(sizeof(mem_chain_v)*batch);
    mem_alnreg_v *local_regvs = malloc(sizeof(mem_alnreg_v)*batch);
    /*
     NEO:
     local_chnnvs saved chains for batch of reads.
        every read has a chnv
        every chains have chains.n chain
            every chain have chain.n seeds
     local_regvs save aln result for batch of reads,
        every read has a regv
     */
    
    //chaining
    chainging_batch(w->opt, w->bwt, w->bns, w->pac, w->seqs+start, w->aux[tid],  batch, local_chnvs);
    
    //extension
    seed_extension_batch(w->opt, w->bwt, w->bns, w->pac, w->seqs+start, w->aux[tid],  batch, local_chnvs, local_regvs);
    
    //post extension
    post_extensiong_batch(w->opt, w->bwt, w->bns, w->pac, w->seqs+start, batch, local_regvs, w->regs+start);
    
    
    //finalize
    for(int i=start, j=0; j<batch; j++,i++)
    {
        mem_chain_v chn = local_chnvs[j];
        for (int i = 0; i < chn.n; ++i) {
            free(chn.a[i].seeds);
        }
        free(chn.a);
    }
    free(local_chnvs);
    free(local_regvs);
}
/*********************************************************/
/*********************************************************/
/*this function is modified by Lingqi Zhang*/
void mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0)
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
    extern void kt_for_batch(int n_threads, int batch_size, void (*func)(void*,int,int), void *data, int n);
    extern void kt_for_batch2(int n_threads, int batch_size, void (*func)(void*,int,int,int), void *data, long n);
	worker_t_mod w;
	mem_pestat_t pes[4];
	double ctime, rtime;
	int i;

    int batch = 1000;
    int batch_size = batch*((opt->flag&MEM_F_PE)?2:1) ;
    
	ctime = cputime(); rtime = realtime();
	global_bns = bns;
	w.regs = malloc(n * sizeof(mem_alnreg_v));
 //   w.chn = malloc(batch_size*opt->n_threads * sizeof(mem_chain_v));
 //   w.ext_val=malloc(batch_size*opt->n_threads * sizeof(qext_t));
    
	w.opt = opt; w.bwt = bwt; w.bns = bns; w.pac = pac;
	w.seqs = seqs; w.n_processed = n_processed;
	w.pes = &pes[0];
	w.aux = malloc(opt->n_threads * sizeof(smem_aux_t));
    
	for (i = 0; i < opt->n_threads; ++i)
		w.aux[i] = smem_aux_init();
    /*
     NEO: this statement need to change to batch mode
     */
    if (bwa_verbose >= 4) printf("=====> Processing %d batchs of read <=====\n", n);
   // kt_for_batch(opt->n_threads, (opt->flag&MEM_F_PE)?2:1, worker_gen_chains, &w, n);
    //fprintf(stderr,"=====> size of swrst_t is %ld   <=====\n", sizeof(swrst_t));
    fprintf(stderr,"=====> Processing %d batchs of read <=====\n", n);
    fprintf(stderr,"=====> Processing %d batchs of read iner <=====\n", batch_size);
    kt_for_batch2(opt->n_threads, batch_size, worker1_batch, &w, n);
    
  //  free(w.ext_val);
#ifdef DEBUG
    printcount();
    reset();
#endif
//    kt_for_batch(opt->n_threads, 10*(opt->flag&MEM_F_PE)?2:1, worker_aln2regs, &w, n);
    
    
    for (i = 0; i < opt->n_threads; ++i)
		smem_aux_destroy(w.aux[i]);
	free(w.aux);
	if (opt->flag&MEM_F_PE) { // infer insert sizes if not provided
		if (pes0) memcpy(pes, pes0, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size distribution as pes0
		else mem_pestat(opt, bns->l_pac, n, w.regs, pes); // otherwise, infer the insert size distribution from data
	}
	kt_for(opt->n_threads, worker2, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // generate alignment
 //   free(w.chn);
//    for(int i=0; i<n; i++)
//    {
//        free(w.regs[i].a);
//    }
    free(w.regs);
	if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);
}
/*********************************************************/
