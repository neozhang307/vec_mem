#include"ksw_batch_simd.h"
#include <stdlib.h>

#include <assert.h>
#include <emmintrin.h>
#include<stdio.h>
#include<string.h>

#include"ksw.h"

void store(swrst_t* data, size_t size, const char* filename)
{
    FILE* output = fopen(filename,"wb+");
    assert(output!=NULL);
    fwrite(&size, sizeof(size_t), 1, output);
    for(int i=0; i<size; i++)
    {
        fwrite(data, sizeof(swrst_t)*size, 1, output);
    }
  
    for(int i=0; i<size; i++)
    {
        int qlen = data[i].sw_seq->qlen;
        
        int rlen = data[i].sw_seq->rlen;
        fwrite(&qlen, sizeof(int), 1, output);
        fwrite(&rlen, sizeof(int), 1, output);
    }
    for(int i=0; i<size; i++)
    {
        
        int qlen = data[i].sw_seq->qlen;
        int rlen = data[i].sw_seq->rlen;
        const uint8_t* q = data[i].sw_seq->query;
        const uint8_t* r = data[i].sw_seq->ref;
        if(qlen!=0)
        {

            fwrite(q, sizeof(uint8_t)*qlen, 1, output);
        }
        if(rlen!=0)
        {
            fwrite(r, sizeof(uint8_t)*rlen, 1, output);
        }

    }
    fclose(output);
}

size_t load(swrst_t** data, const char* filename)
{
    FILE* input = fopen(filename,"rb");
    if(input==NULL)
        return -1;
    size_t size;
    fread(&size, sizeof(size_t), 1, input);
    *data = malloc(sizeof(swrst_t)*size);
    swseq_t* seqs = malloc(sizeof(swseq_t)*size);
    for(int i=0; i<size; i++)
    {
        fread(*data, sizeof(swrst_t)*size, 1, input);
    }
    
    for(int i=0; i<size; i++)
    {
        int qlen;// = data->sw_seq->qlen;
        int rlen;// = data->sw_seq->qlen;
        fread(&qlen, sizeof(int), 1, input);
        fread(&rlen, sizeof(int), 1, input);
        seqs[i].qlen = qlen;
        seqs[i].rlen = rlen;
    //    assert(qlen>=0);
      //  assert(rlen>=0);
    }
    for(int i=0; i<size; i++)
    {
        int qlen = seqs[i].qlen;
        int rlen = seqs[i].rlen;
        uint8_t*q;
        uint8_t*r;
        if(qlen>0)
        {
            q = malloc(sizeof(uint8_t)* qlen);
            fread(q, sizeof(uint8_t)*qlen, 1, input);
        }
        else{
            q = NULL;
        }
        
        if(rlen>0)
        {
            r = malloc(sizeof(uint8_t)* rlen);
            fread(r, sizeof(uint8_t)*rlen, 1, input);
        }
        else{
            r=NULL;
        }
        seqs[i].query = q;
        seqs[i].ref = r;
        (*data)[i].sw_seq=&seqs[i];
    }
    fclose(input);
    return size;
}
static const int g_m = 5;
static int8_t g_mat[5][5];
static int g_o_del;
static int g_e_del;
static int g_o_ins;
static int g_e_ins;
static int g_zdrop;

void finalize_load(swrst_t*data,size_t size)
{
    if(data==NULL)return;
    if(data->sw_seq!=NULL)
    {
        for(int i=0; i<size; i++)
        {
            free((uint8_t*)data[i].sw_seq->query);
            data[i].sw_seq->query=NULL;
            free((uint8_t*)data[i].sw_seq->ref);
            data[i].sw_seq->ref=NULL;
        }
        free(data->sw_seq);
        data->sw_seq = NULL;
    }
    free(data);
}

void load_config()
{
    FILE* confile = fopen("sw.config","rb");
    if(confile==NULL)
        return;
    fread(g_mat, sizeof(int8_t), 25, confile);
    fread(&g_o_del, sizeof(int), 1, confile);
    fread(&g_e_del, sizeof(int), 1, confile);
    fread(&g_o_ins, sizeof(int), 1, confile);
    fread(&g_e_ins, sizeof(int), 1, confile);
    fread(&g_zdrop, sizeof(int), 1, confile);
    fclose(confile);
}

void store_config()
{
    FILE* confile = fopen("sw.config","wb+");
    if(confile==NULL)
        return;
    fwrite(g_mat, sizeof(int8_t), 25, confile);
    fwrite(&g_o_del, sizeof(int), 1, confile);
    fwrite(&g_e_del, sizeof(int), 1, confile);
    fwrite(&g_o_ins, sizeof(int), 1, confile);
    fwrite(&g_e_ins, sizeof(int), 1, confile);
    fwrite(&g_zdrop, sizeof(int), 1, confile);
    fclose(confile);
}

void init(int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int zdrop)
{
    assert(m==5);
    memcpy(g_mat,mat,sizeof(int8_t)*25);
    g_o_del = o_del;
    g_e_del = e_del;
    g_o_ins = o_ins;
    g_e_ins = e_ins;
    g_zdrop = zdrop;
}
int check_config(int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int zdrop)
{
    if(m!=5)
    {
        return -1;
    }
    if(o_del!=g_o_del)
    {
        return -2;
    }
    if(e_del!=g_e_del)
    {
        return -3;
    }
    if(o_ins!=g_o_ins)
    {
        return -5;
    }
    if(e_ins!=g_e_ins)
    {
        return -6;
    }
    if(zdrop!=g_zdrop)
    {
        return -7;
    }
    for(int i=0; i<25; i++)
    {
        if(g_mat[0][i]!=mat[i])
            return -1*i*10;
    }
    return 1;
}
void printcmp(swrst_t* a, swrst_t*b, size_t size);


uint8_t printdif(swrst_t* a, swrst_t*b, size_t size)
{
    uint8_t ret = 0;
    for(int i=0; i<size; i++)
    {
        //printf("process %d of sw\n",i);
    
        swrst_t a_ = a[i];
        swrst_t b_ = b[i];
        if(a_.score!=b_.score)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<0;
            continue;
        }
        if(a_.qle!=b_.qle)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<1;
            continue;
        }if(a_.tle!=b_.tle)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<2;
            continue;
        }
        if(a_.gtle!=b_.gtle)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<3;
            continue;
        }
        if(a_.gscore!=b_.gscore)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<4;
            continue;
        }
        if(a_.max_off!=b_.max_off)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<5;
            continue;
        }
        swseq_t* a_seq = a_.sw_seq;
        swseq_t* b_seq = b_.sw_seq;
        assert(a_seq->qlen==b_seq->qlen);
        for(int i=0; i<a_seq->qlen; i++)
        {
            if(a_seq->query[i]!=b_seq->query[i])
            {
                printcmp(&a_,&b_,1);
                ret|=1<<6;
                break;
            }
        }
        assert(a_seq->rlen==b_seq->rlen);
        for(int i=0; i<a_seq->rlen; i++)
        {
            if(a_seq->ref[i]!=b_seq->ref[i])
            {
                printcmp(&a_,&b_,1);
                ret|=1<<7;
                break;
            }
        }
    }
    return ret;
}

void printcmp(swrst_t* a, swrst_t*b, size_t size)
{
    for(int i=0; i<size; i++)
    {
        swrst_t a_ = a[i];
        swrst_t b_ = b[i];
        
        swseq_t* a_seq = a_.sw_seq;
        swseq_t* b_seq = b_.sw_seq;
        
        printf("A: score: %d, qle: %d, tle:%d, gtle:%d, gscore:%d, max_off:%d qlen:%d rlen:%d\n", a_.score,a_.qle,a_.tle,a_.gtle, a_.gscore, a_.max_off, a_seq->qlen, a_seq->rlen);
        printf("B: score: %d, qle: %d, tle:%d, gtle:%d, gscore:%d, max_off:%d qlen:%d rlen:%d\n", b_.score,b_.qle,b_.tle,b_.gtle, b_.gscore, b_.max_off, b_seq->qlen, a_seq->rlen);
        printf("AQ:\n");
        for(int i=0; i<a_seq->qlen; i++)
        {
            printf("%d",a_seq->query[i]);
        }
        printf("\nBQ:\n");
        for(int i=0; i<b_seq->qlen; i++)
        {
            printf("%d",b_seq->query[i]);
        }
        printf("\nAR:\n");
        for(int i=0; i<a_seq->rlen; i++)
        {
            printf("%d",a_seq->ref[i]);
        }
        printf("\nBR:\n");
        for(int i=0; i<b_seq->rlen; i++)
        {
            printf("%d",b_seq->ref[i]);
        }
        printf("\n");
    }
}

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif


typedef struct {
    int32_t h, e;
} eh_m;



#include"ksort.h"


//ensure the sequence len is ascending
#define sortlen(a, b) ((uint32_t)(a) < (uint32_t)(b))
KSORT_INIT(uint64_t, uint64_t, sortlen)

#define BATCHSIZE 8

typedef struct{
    uint32_t qlen;
    uint32_t rlen;
    uint32_t len;
    uint32_t alined;
    size_t global_batch_id;//index of a batch in whole DB
    uint8_t local_id_y;//idx inside batch
}hash_t;
#define CHECK do{fprintf(stderr,"successfully process to line %d\n",__LINE__);}while(0)
void ksw_extend_batch2(swrst_t* swrts, uint32_t size)
{
    //sort
    assert(size>=0);
    uint64_t* swlen = malloc(sizeof(int64_t)*size);//should record qlen rlen
    for(uint32_t i=0; i<size; ++i)
    {
        uint64_t tval=(uint64_t)i<<32;
        swseq_t* tseq = swrts[i].sw_seq;
        uint32_t qlen = tseq->qlen;
        uint32_t rlen = tseq->rlen;
        uint32_t mx = qlen>rlen?qlen:rlen;
        tval|=mx;
        swlen[i]=tval;
    }
    ks_introsort(uint64_t,size,swlen);
    assert((uint32_t)swlen[0]==0);
    int ptr = 0;
    while((uint32_t)swlen[ptr]==0&&ptr<size) ptr++;
    
    //recompute space size
    //uint8_t batch = BATCHSIZE;//16 | 8
    
    uint32_t resize = size-ptr;
    uint32_t resize_segs = (resize+BATCHSIZE-1)/BATCHSIZE;//skip zero ones
    uint32_t aligned_resize = resize_segs*BATCHSIZE;
    uint64_t *swlen_resized = swlen+ptr;
    
    hash_t* db_hash = malloc(sizeof(hash_t)*aligned_resize);//index related values should be aligned
    memset(db_hash,0,sizeof(hash_t)*aligned_resize);
    
    int global_id_x = 0;
    //break a loop, which is easy to program
    {
        
        int aligned_len = 0;
        int global_batch_id =0;
        int i=0;//i is batch id x;
        for(i=0; i<resize_segs-1; i++)//leave last seg alone
        {
            aligned_len =(uint32_t)swlen_resized[(i+1)*BATCHSIZE-1];//max one is set to use for alignement
            uint64_t * tptr = &swlen_resized[i*BATCHSIZE];
            //last one in a batch
            global_batch_id = global_id_x*BATCHSIZE;
            for(int j=0; j<BATCHSIZE; j++)
            {
                hash_t* cur_rdb = db_hash + i*BATCHSIZE+j;
                cur_rdb->len = (uint32_t)tptr[j];
                cur_rdb->alined = aligned_len;
                cur_rdb->local_id_y = j;
                cur_rdb->global_batch_id = global_batch_id;
            }
            global_id_x+=aligned_len;//point to next position
        }
        aligned_len =(uint32_t)swlen[size-1];//must be the last one
        
        uint64_t * tptr =  &swlen_resized[i*BATCHSIZE];
        int j = 0;
        global_batch_id = global_id_x*BATCHSIZE;
        for(; (i*BATCHSIZE+j)<resize&&j<BATCHSIZE;j++)
        {
            hash_t* cur_rdb = db_hash + i*BATCHSIZE+j;
            cur_rdb->len = (uint32_t)tptr[j];
            cur_rdb->alined = aligned_len;
            cur_rdb->local_id_y = j;
            cur_rdb->global_batch_id = global_batch_id;
        }
        
        for(;j<BATCHSIZE;j++)
        {
            hash_t* cur_rdb = db_hash + i*BATCHSIZE+j;
            cur_rdb->len = 0;
            cur_rdb->alined = aligned_len;
            cur_rdb->local_id_y = j;
            cur_rdb->global_batch_id = global_batch_id;
        }
        global_id_x +=aligned_len;//point to end position
    }
    uint8_t* rdb = malloc(global_id_x*BATCHSIZE);
    memset(rdb,0,sizeof(uint8_t)*global_id_x*BATCHSIZE);
    //copy reference to db
    for(int i=0; i<resize; i++)
    {
        swrst_t *sw = swrts+(swlen_resized[i]>>32);
        swseq_t *seq = sw->sw_seq;
        hash_t* rdb_hash_t = &db_hash[i];
        rdb_hash_t->rlen=seq->rlen;
        
        uint8_t* db_ptr = rdb+rdb_hash_t->global_batch_id+rdb_hash_t->local_id_y*rdb_hash_t->alined;
        memcpy(db_ptr,seq->ref,seq->rlen*sizeof(uint8_t));
    }
    //generate profile
    int8_t* qp_db = malloc(global_id_x*BATCHSIZE*g_m);//should be usingned ,change in the future
    memset(qp_db,(int8_t)-1,sizeof(int8_t)*global_id_x*BATCHSIZE*g_m);
    for(int i=0; i<resize; i++)
    {
        swrst_t *sw = swrts+(swlen_resized[i]>>32);
        swseq_t *seq = sw->sw_seq;
        hash_t* rdb_hash_t = &db_hash[i];
        int qlen =seq->qlen;
        rdb_hash_t->qlen=seq->qlen;
        int8_t* qp = qp_db+g_m*(rdb_hash_t->global_batch_id+rdb_hash_t->local_id_y*rdb_hash_t->alined);
        const uint8_t* query =seq->query;
        for (int k = 0, m = 0; k < g_m; ++k) {
            const int8_t *p = g_mat[k];
            for (int j = 0; j < qlen; ++j) qp[m++] = p[query[j]];
        }
    }
  //  for(int i=ptr; i<size;++i)
    
    int o_del = g_o_del;
    int e_del = g_e_del;
    int o_ins = g_o_ins;
    int e_ins = g_e_ins;
    int zdrop = g_zdrop;
    const int8_t *mat = g_mat[0];
    int m = g_m;
    
    
    
    
    uint64_t* swlen_batch_id = swlen_resized;//resize
    hash_t* db_hash_batch_id = db_hash;//aligned_resize
    
    uint64_t* swlen_nxt_id = swlen_resized;//resize
    hash_t* db_hash_nxt_id = db_hash;//aligned_resize
    
    uint32_t remain = resize;
    uint32_t next_process = BATCHSIZE;
    for(int i=0; i<resize_segs;++i)
    {
        swlen_batch_id = swlen_resized+i* BATCHSIZE;
        db_hash_batch_id = db_hash+i*BATCHSIZE;
        int16_t g_qle[BATCHSIZE];
        int16_t g_tle[BATCHSIZE];
        int16_t g_gtle[BATCHSIZE];
        int16_t g_gscore[BATCHSIZE];
        int16_t g_max_off[BATCHSIZE];
        int16_t g_score[BATCHSIZE];
        int16_t g_h0[BATCHSIZE];// = sw->h0;
        swlen_nxt_id = swlen_batch_id;
        //threoretically wont work
        int j = 0;
        next_process = remain<BATCHSIZE?remain:BATCHSIZE;
        remain-=BATCHSIZE;
        fprintf(stderr,"remaining: %d\n",remain);
        for(j=0; j<BATCHSIZE; j++)
        {
            g_h0[j] = 0;
	}
	for(j=0; j<next_process; j++)
        {
            swrst_t *sw = swrts+((*swlen_nxt_id)>>32);
            g_h0[j] = sw->h0;
            swlen_nxt_id++;
            assert(sw->sw_seq->qlen>0);
            //if(swlen_nxt_id==swlen_end)break;
            //assert(i*BATCHSIZE+j<resize);
        }
        
        
       // swlen_nxt_id = swlen_batch_id;
        db_hash_nxt_id = db_hash_batch_id;
        for(int j=0; j<next_process; j++)
            //process 16 query at a time for uint8; process 8 query at a time for uint16 query
            {
                
                int qlen = db_hash_nxt_id->qlen;
                int tlen = db_hash_nxt_id->rlen;//seq->rlen;
                size_t seed_global_id = db_hash_nxt_id->global_batch_id+db_hash_nxt_id->local_id_y*db_hash_nxt_id->alined;
                const uint8_t *target =  rdb+seed_global_id;//seq->ref;
                int8_t *qp = qp_db+g_m*(seed_global_id);//malloc(qlen * m);
                int ali_len = db_hash_nxt_id->alined;
                
                eh_m *eh; // score array
                 // query profile
                int16_t i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end;
                int16_t max, max_i, max_j, max_ie, gscore, max_off;
                int16_t h0=g_h0[j];
                assert(h0 >= 0);
                // allocate memory
                eh = calloc(qlen + 1, 8);
                
                // fill the first row
                eh[0].h = h0; eh[1].h = h0 > oe_ins? h0 - oe_ins : 0;
                for (j = 2; j <= qlen && eh[j-1].h > e_ins; ++j)
                    eh[j].h = eh[j-1].h - e_ins;
                // adjust $w if it is too large
                k = m * m;
                for (i = 0, max = 0; i < k; ++i) // get the max score
                    max = max > mat[i]? max : mat[i];
                // DP loop
                max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
                max_off = 0;
                beg = 0, end = qlen;
               
                //MAIN SW
                //for (i = 0; LIKELY(i < tlen); ++i) {
                for (i = 0; LIKELY(i < ali_len); ++i) {
                    int t, f = 0, h1, m = 0, mj = -1;
                    int8_t *q = &qp[target[i] * qlen];
                    // apply the band and the constraint (if provided)
                    //        if (beg < i - w) beg = i - w;
                    //        if (end > i + w + 1) end = i + w + 1;
                    //        if (end > qlen) end = qlen;
                    // compute the first column
                    if (beg == 0) {
                        h1 = h0 - (o_del + e_del * (i + 1));
                        if (h1 < 0) h1 = 0;
                    } else h1 = 0;
                    //processing a row
                    for (j = beg; LIKELY(j < end); ++j) {
                        // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
                        // Similar to SSE2-SW, cells are computed in the following order:
                        //   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
                        //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
                        //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
                        eh_m *p = &eh[j];
                        int h, M = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
                        p->h = h1;          // set H(i,j-1) for the next row
                        M = M? M + q[j] : 0;// separating H and M to disallow a cigar like "100M3I3D20M"
                        h = M > e? M : e;   // e and f are guaranteed to be non-negative, so h>=0 even if M<0
                        h = h > f? h : f;
                        h1 = h;             // save H(i,j) to h1 for the next column
                        mj = m > h? mj : j; // record the position where max score is achieved
                        m = m > h? m : h;   // m is stored at eh[mj+1]
                        t = M - oe_del;
                        t = t > 0? t : 0;
                        e -= e_del;
                        e = e > t? e : t;   // computed E(i+1,j)
                        p->e = e;           // save E(i+1,j) for the next row
                        t = M - oe_ins;
                        t = t > 0? t : 0;
                        f -= e_ins;
                        f = f > t? f : t;   // computed F(i,j+1)
                    }
                    h1=i<tlen?h1:0;
                    eh[end].h = h1; eh[end].e = 0;
                    if (j == qlen) {
                        max_ie = gscore > h1? max_ie : i;
                        gscore = gscore > h1? gscore : h1;
                    }
                    if (m == 0) break; //theoretically not important , can be change to bach
                    if (m > max) {
                        max = m, max_i = i, max_j = mj;
                        max_off = max_off > abs(mj - i)? max_off : abs(mj - i);
                    } else if (zdrop > 0) {
                        if (i - max_i > mj - max_j) {
                            if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break;
                        } else {
                            if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break;
                        }
                    }
                    // update beg and end for the next round
//                    for (j = beg; LIKELY(j < end) && eh[j].h == 0 && eh[j].e == 0; ++j);
//                    beg = j;
//                    for (j = end; LIKELY(j >= beg) && eh[j].h == 0 && eh[j].e == 0; --j);
//                    end = j + 2 < qlen? j + 2 : qlen;
                    //beg = 0; end = qlen; // uncomment this line for debugging
                }
                
                //finalize
                free(eh); //free(qp);

                //post process
                //result
            
                g_qle[j] = max_j+1;
                g_tle[j] = max_i+1;
                g_gtle[j] = max_ie+1;
                g_gscore[j] = gscore;
                g_max_off[j] = max_off;
                g_score[j] = max;
                //swlen_nxt_id++;
                db_hash_nxt_id++;
                
            }
        
        swlen_nxt_id = swlen_batch_id;
        for(j=0; j<next_process; j++)
        {
            swrst_t *sw = swrts+((*swlen_nxt_id)>>32);
            sw->qle = g_qle[j];
            sw->tle = g_tle[j];
            sw->gtle = g_gtle[j];
            sw->gscore = g_gscore[j];
            sw->max_off = g_max_off[j];
            sw->score = g_score[j];
            swlen_nxt_id++;
          //  if(swlen_nxt_id==swlen_end)break;
        }
        
//        swlen_batch_id+=BATCHSIZE;
//        db_hash_batch_id+=BATCHSIZE;
    }
    
    
    free(qp_db);
    free(rdb);
    free(swlen);
    free(db_hash);
}



void ksw_extend_batch(swrst_t* swrts, size_t size)
{
    //sort
    
    for(int i=0; i<size;++i)
    {
//            swrst_t *sw = swrts+i;
//            swseq_t *seq = sw->sw_seq;
//            if(seq->qlen!=0)
//            {
//                sw->score = ksw_extend2_mod(seq->qlen, seq->query, seq->rlen,seq->ref, 5, g_mat[0], g_o_del, g_e_del, g_o_ins, g_e_ins, g_zdrop, sw->h0, &sw->qle, &sw->tle, &sw->gtle, &sw->gscore, &sw->max_off);
//            }
        swrst_t *sw = swrts+i;
        swseq_t *seq = sw->sw_seq;
        if(seq->qlen!=0)
        {
            
            int qlen = seq->qlen;
            int tlen = seq->rlen;
            const uint8_t *query = seq->query;
            const uint8_t *target = seq->ref;
            int o_del = g_o_del;
            int e_del = g_e_del;
            int o_ins = g_o_ins;
            int e_ins = g_e_ins;
            int zdrop = g_zdrop;
            int h0 = sw->h0;
            const int8_t *mat = g_mat[0];
            int m = g_m;
            //process 16 query at a time for uint8; process 8 query at a time for uint16 query
            {
                eh_m *eh; // score array
                int8_t *qp; // query profile
                int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end, max, max_i, max_j, max_ie, gscore, max_off;
                assert(h0 > 0);
                // allocate memory
                qp = malloc(qlen * m);
                eh = calloc(qlen + 1, 8);
                // generate the query profile
                for (k = i = 0; k < m; ++k) {
                    const int8_t *p = &mat[k * m];
                    for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
                }
                
                
                // fill the first row
                eh[0].h = h0; eh[1].h = h0 > oe_ins? h0 - oe_ins : 0;
                for (j = 2; j <= qlen && eh[j-1].h > e_ins; ++j)
                    eh[j].h = eh[j-1].h - e_ins;
                // adjust $w if it is too large
                k = m * m;
                for (i = 0, max = 0; i < k; ++i) // get the max score
                    max = max > mat[i]? max : mat[i];
                // DP loop
                max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
                max_off = 0;
                beg = 0, end = qlen;
                
                //MAIN SW
                for (i = 0; LIKELY(i < tlen); ++i) {
                    int t, f = 0, h1, m = 0, mj = -1;
                    int8_t *q = &qp[target[i] * qlen];
                    // apply the band and the constraint (if provided)
                    //        if (beg < i - w) beg = i - w;
                    //        if (end > i + w + 1) end = i + w + 1;
                    //        if (end > qlen) end = qlen;
                    // compute the first column
                    if (beg == 0) {
                        h1 = h0 - (o_del + e_del * (i + 1));
                        if (h1 < 0) h1 = 0;
                    } else h1 = 0;
                    //processing a row
                    for (j = beg; LIKELY(j < end); ++j) {
                        // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
                        // Similar to SSE2-SW, cells are computed in the following order:
                        //   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
                        //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
                        //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
                        eh_m *p = &eh[j];
                        int h, M = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
                        p->h = h1;          // set H(i,j-1) for the next row
                        M = M? M + q[j] : 0;// separating H and M to disallow a cigar like "100M3I3D20M"
                        h = M > e? M : e;   // e and f are guaranteed to be non-negative, so h>=0 even if M<0
                        h = h > f? h : f;
                        h1 = h;             // save H(i,j) to h1 for the next column
                        mj = m > h? mj : j; // record the position where max score is achieved
                        m = m > h? m : h;   // m is stored at eh[mj+1]
                        t = M - oe_del;
                        t = t > 0? t : 0;
                        e -= e_del;
                        e = e > t? e : t;   // computed E(i+1,j)
                        p->e = e;           // save E(i+1,j) for the next row
                        t = M - oe_ins;
                        t = t > 0? t : 0;
                        f -= e_ins;
                        f = f > t? f : t;   // computed F(i,j+1)
                    }
                    eh[end].h = h1; eh[end].e = 0;
                    if (j == qlen) {
                        max_ie = gscore > h1? max_ie : i;
                        gscore = gscore > h1? gscore : h1;
                    }
                    if (m == 0) break;
                    if (m > max) {
                        max = m, max_i = i, max_j = mj;
                        max_off = max_off > abs(mj - i)? max_off : abs(mj - i);
                    } else if (zdrop > 0) {
                        if (i - max_i > mj - max_j) {
                            if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break;
                        } else {
                            if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break;
                        }
                    }
                    // update beg and end for the next round
         //           for (j = beg; LIKELY(j < end) && eh[j].h == 0 && eh[j].e == 0; ++j);
          //          beg = j;
           //         for (j = end; LIKELY(j >= beg) && eh[j].h == 0 && eh[j].e == 0; --j);
            //        end = j + 2 < qlen? j + 2 : qlen;
                    beg = 0; end = qlen; // uncomment this line for debugging
                }
                
                //finalize
                free(eh); free(qp);
                
                //post process
                //result
                sw->qle = max_j+1;
                sw->tle = max_i+1;
                sw->gtle = max_ie+1;
                sw->gscore = gscore;
                sw->max_off = max_off;
                sw->score = max;
            }
        }
    }
}
#ifdef SWBATCHDB
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>

double cputime()
{
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

int main()
{
    load_config();
    swrst_t* nsrt;
    size_t nread = load(&nsrt,"sw_start_8000_0_2000.bin");
    
    //time
    double ctime = cputime();
    double rtime = realtime();
    
    ksw_extend_batch2(nsrt, nread);
    
    //time
    fprintf(stderr, "[M::%s] Processed %ld reads in %.3f CPU sec, %.3f real sec\n", __func__, nread, cputime() - ctime, realtime() - rtime);

    swrst_t* rsrt;
    size_t rread = load(&rsrt,"sw_end_8000_0_2000.bin");
    assert(rread==nread);
    //printf("the result is %d\n", cmp(nsrt,rsrt,nread));
    uint8_t check = printdif(nsrt,rsrt,nread);
    if(check==0)
    {
        fprintf(stderr,"the check result is correct\n");
        fprintf(stderr,"which is 0x%02x\n",check);
    }
    else
    {
        fprintf(stderr,"the check result is incorrect\n");
        fprintf(stderr,"which is 0x%02x\n",check);
    }
    finalize_load(nsrt,nread);
    finalize_load(rsrt,rread);
    return 0;
}
#endif
