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
   #define print(i) do{\
            \
            {\
             fprintf(stderr,"i:%d ",i);   \
                \
            }\
    }while(0)
    for(int i=0; i<size; i++)
    {
        //printf("process %d of sw\n",i);
    
        swrst_t a_ = a[i];
        swrst_t b_ = b[i];
        if(a_.score!=b_.score)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<0;
            print(i);
            continue;
        }
        if(a_.qle!=b_.qle)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<1;
            print(i);
            continue;
        }if(a_.tle!=b_.tle)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<2;
            print(i);
            continue;
        }
        if(a_.gtle!=b_.gtle)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<3;
            print(i);
            continue;
        }
        if(a_.gscore!=b_.gscore)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<4;
            print(i);
            continue;
        }
        if(a_.max_off!=b_.max_off)
        {
            printcmp(&a_,&b_,1);
            ret|=1<<5;
            print(i);
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
                print(i);
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
                print(i);
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
#include"ktranspose.h"


//ensure the sequence len is ascending
#define sortlen(a, b) ((uint32_t)(a) < (uint32_t)(b))
KSORT_INIT(uint64_t, uint64_t, sortlen)

#define BATCHSIZE 16
#define PROCESSBATCH 8
typedef struct{
    uint32_t qlen;
    uint32_t rlen;
    uint32_t len;
    uint32_t alined;
    size_t global_batch_id;//index of a batch in whole DB
    uint8_t local_id_y;//idx inside batch
}hash_t;
#define CHECK do{fprintf(stderr,"successfully process to line %d\n",__LINE__);}while(0)
static int one=0;
void batch_sw_core2(hash_t* db_hash_batch_id,
                   uint8_t* rdb_rev,
                   int16_t* qp_db,
                   int16_t g_h0[BATCHSIZE],//input
                   int o_del,
                   int e_del,
                   int o_ins,
                   int e_ins,
                   int zdrop,
                   
                   int16_t g_qle[BATCHSIZE],//result
                   int16_t g_tle[BATCHSIZE],
                   int16_t g_gtle[BATCHSIZE],
                   int16_t g_gscore[BATCHSIZE],
                   int16_t g_max_off[BATCHSIZE],
                   int16_t g_score[BATCHSIZE]
                   
                   )
{
    if(one!=1){
        fprintf(stderr, "processing new SW core\n");
        one=1;
    }
#define __max_8(ret, xx) do { \
(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
(ret) = _mm_extract_epi16((xx), 0); \
} while (0)
#define __min_8(ret, xx) do { \
(xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 8)); \
(xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 4)); \
(xx) = _mm_min_epi16((xx), _mm_srli_si128((xx), 2)); \
(ret) = _mm_extract_epi16((xx), 0); \
} while (0)
    __m128i v_zero, v_oe_del, v_e_del, v_oe_ins, v_e_ins;
    
    v_zero = _mm_setzero_pd();
    v_oe_del = _mm_set1_epi16(o_del + e_del);
    v_e_del = _mm_set1_epi16(e_del);
    v_oe_ins = _mm_set1_epi16(o_ins + e_ins);
    v_e_ins = _mm_set1_epi16(e_ins);
    
    int64_t  oe_del = o_del + e_del, oe_ins = o_ins + e_ins;
    
    uint16_t p_tlen[8],p_qlen[8];
    __m128i v_tlen, v_qlen;
    __m128i v_beg, v_end, v_align_end;
    hash_t*  db_hash_nxt_id = db_hash_batch_id;
    
    __m128i v_max, v_max_i, v_max_j, v_max_ie, v_gscore;
    __m128i v_max_off;
    __m128i v_H0, v_H1;
  //  max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
  //  max_off = 0;
    for(int grid_process_batch_idx=0; grid_process_batch_idx<BATCHSIZE/PROCESSBATCH;grid_process_batch_idx++)
    {
        //process 8 query at a time for int16_t
        
        int ali_len = db_hash_nxt_id->alined;
        for(int process_batch_id = 0; process_batch_id<PROCESSBATCH; process_batch_id++)
        {
            
            p_qlen[process_batch_id]=db_hash_nxt_id[process_batch_id].qlen;
            p_tlen[process_batch_id]=db_hash_nxt_id[process_batch_id].rlen;
        }
        v_qlen = _mm_load_si128((__m128i*)p_qlen);//theoretically not necessary
        v_tlen = _mm_load_si128((__m128i*)p_tlen);
        v_beg=_mm_setzero_pd();
        v_end = _mm_load_si128((__m128i*)p_qlen);
        v_align_end = _mm_set1_epi16(ali_len);
        
        __m128i v_qp[PROCESSBATCH];//8X8
        int16_t *qp = qp_db+g_m*(db_hash_nxt_id[0].global_batch_id+(0 + grid_process_batch_idx*PROCESSBATCH)*ali_len);
//        int16_t *qp1 = qp_db+g_m*(db_hash_nxt_id[1].global_batch_id+(1 + grid_process_batch_idx*PROCESSBATCH)*ali_len);
//        int16_t *qp2 = qp_db+g_m*(db_hash_nxt_id[2].global_batch_id+(2 + grid_process_batch_idx*PROCESSBATCH)*ali_len);
//        int16_t *qp3 = qp_db+g_m*(db_hash_nxt_id[3].global_batch_id+(3 + grid_process_batch_idx*PROCESSBATCH)*ali_len);
//        int16_t *qp4 = qp_db+g_m*(db_hash_nxt_id[4].global_batch_id+(4 + grid_process_batch_idx*PROCESSBATCH)*ali_len);
//        int16_t *qp5 = qp_db+g_m*(db_hash_nxt_id[5].global_batch_id+(5 + grid_process_batch_idx*PROCESSBATCH)*ali_len);
//        int16_t *qp6 = qp_db+g_m*(db_hash_nxt_id[6].global_batch_id+(6 + grid_process_batch_idx*PROCESSBATCH)*ali_len);
//        int16_t *qp7 = qp_db+g_m*(db_hash_nxt_id[7].global_batch_id+(7 + grid_process_batch_idx*PROCESSBATCH)*ali_len);
        const uint8_t *target_rev_batch =  rdb_rev+db_hash_nxt_id[0].global_batch_id;
        
        __m128i *H,*E;
        //q0 q1 q2 q3 q4 q5 q6 q7
        //q0 q1 q2 q3 q4 q5 q6 q7
        //...
        // count = ali_len+1
        H = (__m128i*)malloc(sizeof(__m128i)*(ali_len + 1));
        memset(H,0,sizeof(__m128i)*(ali_len + 1));
        E = (__m128i*)malloc(sizeof(__m128i)*(ali_len + 1));
        memset(E,0,sizeof(__m128i)*(ali_len + 1));
        v_H0 = _mm_load_si128(((__m128i*)g_h0)+grid_process_batch_idx);
        H[0]=v_H0;
        H[1] = _mm_subs_epu16(H[0],v_oe_ins);
        for(int itr_j = 2; itr_j<=ali_len; itr_j++)
        {
            H[itr_j]=_mm_subs_epu16(H[itr_j-1],v_e_ins);
        }
        
        v_max = v_H0;
        v_gscore = _mm_set1_epi16(-1);
        v_max_ie = _mm_set1_epi16(-1);
        v_max_i=_mm_set1_epi16(-1);
        v_max_j=_mm_set1_epi16(-1);
        v_max_off = _mm_setzero_pd();
        int i_end;
        __max_8(i_end, v_tlen);
        for (int i = 0; LIKELY(i < i_end); ++i) {//tlen
//            int t, f = 0, h1, m = 0, mj = -1;
//            int h_l=0, m_l=0, mj_l=0;
            __m128i v_t, v_f, v_m, v_mj, v_h_l, v_m_l, v_mj_l;
            
            //init qprofile
            const uint8_t *nxt_target = target_rev_batch+i*BATCHSIZE+grid_process_batch_idx*PROCESSBATCH;
            uint8_t v_target[8];
            memcpy(v_target,nxt_target,8);
            
            
            v_f = _mm_setzero_pd();
            v_m = _mm_setzero_pd();
            v_mj = _mm_set1_epi16(-1);
            v_h_l = _mm_setzero_pd();
            v_m_l = _mm_setzero_pd();
            v_mj_l = _mm_setzero_pd();
            //            if (beg == 0) {
            //                h1 = h0 - (o_del + e_del * (i + 1));
            //                if (h1 < 0) h1 = 0;
            //            } else h1 = 0;
            __m128i cmp =_mm_cmpeq_epi16(v_beg, v_zero);
            
            __m128i tmp =_mm_set1_epi16(o_del+e_del*(i+1));
            v_H1 = _mm_subs_epu16(v_H0,tmp);//if (beg == 0) {h1 = h0 - (o_del + e_del * (i + 1));
            v_H1 = _mm_and_si128(v_H1, cmp);//} else h1 = 0;
//            for (j = beg; LIKELY(j < align_end); ++j) {
            int j_beg, j_end;
            __min_8(j_beg,v_beg);
            __max_8(j_end,v_end);
            for(int j=j_beg; j<j_end; j++)//end
            {
               
                
            }
            
            break;
        }
        for(int process_batch_id = 0; process_batch_id<PROCESSBATCH; process_batch_id++)
        {
            int batch_idx = process_batch_id + grid_process_batch_idx*PROCESSBATCH;
            int qlen = db_hash_nxt_id->qlen;
            int tlen = db_hash_nxt_id->rlen;//seq->rlen;
           // CHECK here assert(p_qlen[process_batch_id] == qlen);
            /************************/
            size_t batch_global_id =  db_hash_nxt_id->global_batch_id;
            // int alignLen = db_hash_nxt_id->alined;
            
            size_t seed_global_id = db_hash_nxt_id->global_batch_id+batch_idx*db_hash_nxt_id->alined;
            /***********************/
            //const uint8_t *target_batch =  rdb+batch_global_id;//seq->ref; //BATCHSIZE * alignLen;
            const uint8_t *target_rev_batch =  rdb_rev+batch_global_id;//seq->ref; //alignLen * BATCHSIZE
            /***********************/
            int16_t *qp = qp_db+g_m*(seed_global_id);//malloc(qlen * m);
            int ali_len = db_hash_nxt_id->alined;
            
            eh_m *eh; // score array
            // query profile
            int16_t i, j;
            int16_t beg, end, align_end;
            int16_t max, max_i, max_j, max_ie, gscore, max_off;
            int16_t h0=g_h0[batch_idx];
            assert(h0 >= 0);
            // allocate memory
            align_end = ali_len;
            beg = 0, end = qlen;//every seqs in a batch should have same qlen
            
            eh = calloc(align_end + 1, 8);
            
            // fill the first row
            eh[0].h = h0; eh[1].h = h0 > oe_ins? h0 - oe_ins : 0;
            for (j = 2; j <= align_end && eh[j-1].h > e_ins; ++j)
                eh[j].h = eh[j-1].h - e_ins;
#ifdef DEBUG
            for(int i=0; i<ali_len; i++)
            {
                //fprintf(stderr,"%d:%d \t",eh[i].h,((int16_t*)H)[process_batch_id+i*PROCESSBATCH]);
                assert(eh[i].h==((int16_t*)H)[process_batch_id+i*PROCESSBATCH]);
            }
#endif
            // fprintf(stderr,"\n");
            
            // adjust $w if it is too large
            // DP loop
            max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
            max_off = 0;
            //MAIN SW
            for (i = 0; LIKELY(i < ali_len); ++i) {//tlen
                int t, f = 0, h1, m = 0, mj = -1;
                int h_l=0, m_l=0, mj_l=0;
                uint8_t nxt_target = target_rev_batch[i*BATCHSIZE+batch_idx];
                int16_t *q = &qp[nxt_target * ali_len];
                ////////////////////////
                // apply the band and the constraint (if provided)
                //        if (beg < i - w) beg = i - w;
                //        if (end > i + w + 1) end = i + w + 1;
                //        if (end > qlen) end = qlen;
                // compute the first column
                if (beg == 0) {
                    h1 = h0 - (o_del + e_del * (i + 1));
                    if (h1 < 0) h1 = 0;
                } else h1 = 0;
#ifdef DEBUG
                if(i==0)
                {
                    uint16_t a[8];
                    _mm_store_si128((__m128i*)a,v_H1);
                    assert(h1==a[process_batch_id]);
                }
#endif
                //processing a row
                for (j = beg; LIKELY(j < align_end); ++j) {//end
                    // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
                    // Similar to SSE2-SW, cells are computed in the following order:
                    //   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
                    //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
                    //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
                    eh_m *p = &eh[j];
                    int h, M = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
                    p->h = h1;          // set H(i,j-1) for the next row
                    M = M? M + q[j] : 0;// separating H and M to disallow a cigar like "100M3I3D20M"
                    /////////
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
                    //record last h1, m, mj
                    h_l=j<end?h1:h_l;
                    m_l=j<end?m:m_l;
                    mj_l=j<end?mj:mj_l;
                }
                j=j<end?j:end;
                h1=h_l;
                m=m_l;
                mj=mj_l;
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
            
            g_qle[batch_idx] = max_j+1;
            g_tle[batch_idx] = max_i+1;
            g_gtle[batch_idx] = max_ie+1;
            g_gscore[batch_idx] = gscore;
            g_max_off[batch_idx] = max_off;
            g_score[batch_idx] = max;
            db_hash_nxt_id++;
            
        }
    }
}


void batch_sw_core(hash_t* db_hash_batch_id,
                   uint8_t* rdb_rev,
                   int16_t* qp_db,
                   int16_t g_h0[BATCHSIZE],//input
                   int o_del,
                   int e_del,
                   int o_ins,
                   int e_ins,
                   int zdrop,
                   
                   int16_t g_qle[BATCHSIZE],//result
                   int16_t g_tle[BATCHSIZE],
                   int16_t g_gtle[BATCHSIZE],
                   int16_t g_gscore[BATCHSIZE],
                   int16_t g_max_off[BATCHSIZE],
                   int16_t g_score[BATCHSIZE]
                   
                   )
{
    
    hash_t*  db_hash_nxt_id = db_hash_batch_id;
    int16_t oe_del = o_del + e_del, oe_ins = o_ins + e_ins;
    int align_end = db_hash_nxt_id->alined;
    size_t batch_global_id =  db_hash_nxt_id->global_batch_id;

    
    int16_t* qp_buff = malloc(sizeof(int16_t)*PROCESSBATCH*align_end);
    memset(qp_buff,0,sizeof(int16_t)*PROCESSBATCH*align_end);
    
    for(int grid_process_batch_idx=0; grid_process_batch_idx<BATCHSIZE/PROCESSBATCH;grid_process_batch_idx++)
    {
        const uint8_t *target_rev_batch =  rdb_rev+batch_global_id;
        const int16_t *qp_batch = qp_db +g_m*batch_global_id;
        const int16_t *qp_batch_nxt = qp_batch + g_m*grid_process_batch_idx*8*align_end;
        //process 8 query at a time for int16_t
        uint16_t qlens[8];
        uint16_t maxqlen=0;
        uint16_t tlens[8];
        uint16_t maxtlen=0;
        
        
        
        eh_m **ehs=malloc(sizeof(eh_m*)*(align_end + 1));
        for(int i=0; i<align_end+1;i++)
        {
            ehs[i]=calloc(sizeof(eh_m)*8,1);
        }
        int16_t begs[8], ends[8];
        int16_t maxs[8], max_is[8], max_js[8], max_ies[8], gscores[9], max_offs[8];
        int16_t h0s[8];
        for(int process_batch_id=0; process_batch_id<8; process_batch_id++)
        {
            h0s[process_batch_id]=g_h0[process_batch_id+grid_process_batch_idx*8];
            assert( h0s[process_batch_id] >= 0);
        }
        for(int process_batch_id=0; process_batch_id<8; process_batch_id++)
        {
            ehs[0][process_batch_id].h=h0s[process_batch_id];
            ehs[1][process_batch_id].h=h0s[process_batch_id]>oe_ins?h0s[process_batch_id]-oe_ins:0;
            for (int j = 2; j <= align_end && ehs[j-1][process_batch_id].h > e_ins; ++j)
                ehs[j][process_batch_id].h = ehs[j-1][process_batch_id].h - e_ins;
        }
        for(int process_batch_id=0; process_batch_id<8; process_batch_id++)
        {
            hash_t *tmphash =db_hash_nxt_id+process_batch_id+grid_process_batch_idx*8;
            
            qlens[process_batch_id]=tmphash->qlen;
            maxqlen=maxqlen>tmphash->qlen?maxqlen:tmphash->qlen;
            
            tlens[process_batch_id]=tmphash->rlen;
            maxtlen=maxtlen>tmphash->rlen?maxtlen:tmphash->rlen;
        }
        for(int process_batch_id=0; process_batch_id<8; process_batch_id++)
        {
            int batch_idx = process_batch_id + grid_process_batch_idx*8;
            h0s[process_batch_id]=g_h0[batch_idx];
            maxs[process_batch_id]=h0s[process_batch_id];
            max_is[process_batch_id]=-1;
            max_js[process_batch_id]=-1;
            max_ies[process_batch_id]=-1;
            gscores[process_batch_id]=-1;
            max_offs[process_batch_id]=0;
        }
      //  int16_t begs[8];
        //int16_t neds[8];
        int16_t ts[8], fs[8] , h1s[8], ms[8], mjs[8];
        int16_t h_ls[8], m_ls[8], mj_ls[8];
        int16_t min_beg, max_end;
        for(int process_batch_id = 0; process_batch_id<8; process_batch_id++)
        {
            begs[process_batch_id]=0;
            ends[process_batch_id]=qlens[process_batch_id];

        }
        


            //MAIN SW
        uint8_t break_flag = 0;
        for (int16_t i = 0; LIKELY(i < maxtlen) && break_flag==0; ++i) {
            min_beg = begs[0];
            max_end = ends[0];
            for(int i=1; i<8; i++){
                int16_t tbeg = begs[i];
                int16_t tend = ends[i];
                min_beg = min_beg<tbeg?min_beg:tbeg;
                max_end = max_end>tend?max_end:tend;
            }
            
            //int16_t* qp_buff = malloc(sizeof(int16_t)*PROCESSBATCH*align_end);
            //int16_t* nxt_qp_buff = qp_buff;
            uint8_t t_targets[8];
            memcpy(t_targets,target_rev_batch+i*BATCHSIZE + grid_process_batch_idx*8,8*sizeof(uint8_t));
            int16_t* qp_buff_nxt = qp_buff;
            for(int process_batch_id=0; process_batch_id<PROCESSBATCH; process_batch_id++)
            {
                int qp_ptr = process_batch_id*g_m*align_end+t_targets[process_batch_id] * align_end;
                
                memcpy(qp_buff_nxt,qp_batch_nxt+qp_ptr,align_end*sizeof(int16_t));
                qp_buff_nxt+=align_end;
            }

            for(int process_batch_id = 0; process_batch_id<8; process_batch_id++)
            {
                /***********************/
                
                int16_t j;

                fs[process_batch_id]=0;
                ms[process_batch_id]=0;
                mjs[process_batch_id]=-1;
                h_ls[process_batch_id]=0;
                m_ls[process_batch_id]=0;
                mj_ls[process_batch_id]=0;
                
                //int16_t m = ms[process_batch_id];
                //int16_t mj = mjs[process_batch_id];
                //int16_t h_l=h_ls[process_batch_id];
                //int16_t m_l=m_ls[process_batch_id];
                //int16_t mj_l=mj_ls[process_batch_id];
                /***********************/

                const  int16_t *q = qp_buff+align_end*process_batch_id;
                // compute the first column
                if ( min_beg == 0) {
                    h1s[process_batch_id] = h0s[process_batch_id] - (o_del + e_del * (i + 1));
                    if (h1s[process_batch_id] < 0) h1s[process_batch_id] = 0;
                } else h1s[process_batch_id] = 0;
                //processing a row
                for (j =  min_beg; LIKELY(j < max_end); ++j) {
                    // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
                    // Similar to SSE2-SW, cells are computed in the following order:
                    //   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
                    //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
                    //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape

                    int16_t h;
                    int16_t M = ehs[j][process_batch_id].h;
                    int16_t e = ehs[j][process_batch_id].e; // get H(i-1,j-1) and E(i-1,j)
                    
                    ehs[j][process_batch_id].h = h1s[process_batch_id];          // set H(i,j-1) for the next row
                    
                    M = M? M + q[j] : 0;// separating H and M to disallow a cigar like "100M3I3D20M"
                    
                    h = M > e? M : e;   // e and f are guaranteed to be non-negative, so h>=0 even if M<0
                    h = h > fs[process_batch_id]? h : fs[process_batch_id];
                    h1s[process_batch_id] = h;             // save H(i,j) to h1 for the next column
                    mjs[process_batch_id] = ms[process_batch_id] > h?  mjs[process_batch_id] : j; // record the position where max score is achieved
                    ms[process_batch_id] = ms[process_batch_id] > h? ms[process_batch_id] : h;   // m is stored at eh[mj+1]
                    ts[process_batch_id] = M - oe_del;
                    ts[process_batch_id] = ts[process_batch_id] > 0? ts[process_batch_id] : 0;
                    e -= e_del;
                    e = e > ts[process_batch_id]? e : ts[process_batch_id];   // computed E(i+1,j)
                    ehs[j][process_batch_id].e = e;           // save E(i+1,j) for the next row
                    ts[process_batch_id] = M - oe_ins;
                    ts[process_batch_id] = ts[process_batch_id] > 0? ts[process_batch_id] : 0;
                    fs[process_batch_id] -= e_ins;
                    fs[process_batch_id] = fs[process_batch_id] > ts[process_batch_id]? fs[process_batch_id] : ts[process_batch_id];   // computed F(i,j+1)
                    //record last h1, m, mj
                    h_ls[process_batch_id]=j<ends[process_batch_id]?h1s[process_batch_id]:h_ls[process_batch_id];
                    m_ls[process_batch_id]=j<ends[process_batch_id]?ms[process_batch_id]: m_ls[process_batch_id];
                    mj_ls[process_batch_id]=j<ends[process_batch_id]? mjs[process_batch_id]:mj_ls[process_batch_id];
                }
                j=j<ends[process_batch_id]?j:ends[process_batch_id];
                h1s[process_batch_id]=h_ls[process_batch_id];
                ms[process_batch_id]= m_ls[process_batch_id];
                mjs[process_batch_id]=mj_ls[process_batch_id];
                h1s[process_batch_id]=i<tlens[process_batch_id]?h1s[process_batch_id]:0;
                ehs[ends[process_batch_id]][process_batch_id].h = h1s[process_batch_id]; ehs[ends[process_batch_id]][process_batch_id].e = 0;
                if (j == qlens[process_batch_id]) {
                    max_ies[process_batch_id] = gscores[process_batch_id] > h1s[process_batch_id]? max_ies[process_batch_id] : i;
                    gscores[process_batch_id] = gscores[process_batch_id] > h1s[process_batch_id]? gscores[process_batch_id] : h1s[process_batch_id];
                }
   //             if (ms[process_batch_id] == 0) break; //theoretically not important , can be change to bach
//                if (ms[process_batch_id] > maxs[process_batch_id]) {
//                    maxs[process_batch_id] = ms[process_batch_id], max_is[process_batch_id] = i, max_js[process_batch_id] =  mjs[process_batch_id];
//                    max_offs[process_batch_id] = max_offs[process_batch_id] > abs( mjs[process_batch_id] - i)? max_offs[process_batch_id] : abs( mjs[process_batch_id] - i);
//                } else if (zdrop > 0) {
//                    if (i - max_is[process_batch_id] >  mjs[process_batch_id] - max_js[process_batch_id]) {
//                        if (maxs[process_batch_id] - ms[process_batch_id] - ((i - max_is[process_batch_id]) - ( mjs[process_batch_id] - max_js[process_batch_id])) * e_del > zdrop) break;
//                    } else {
//                        if (maxs[process_batch_id] - ms[process_batch_id]- (( mjs[process_batch_id] - max_js[process_batch_id]) - (i - max_is[process_batch_id])) * e_ins > zdrop) break;
//                    }
//                }
            }
            //if the search should terminated earlier?
            uint8_t flag = 0;
            for(int process_batch_id = 0; process_batch_id<8; process_batch_id++)
            {
                // if (ms[process_batch_id] == 0) break;
                if(ms[process_batch_id]!=0)
                    flag=1;
                if (ms[process_batch_id] > maxs[process_batch_id]) {
                    maxs[process_batch_id] = ms[process_batch_id], max_is[process_batch_id] = i, max_js[process_batch_id] =  mjs[process_batch_id];
                    max_offs[process_batch_id] = max_offs[process_batch_id] > abs( mjs[process_batch_id] - i)? max_offs[process_batch_id] : abs( mjs[process_batch_id] - i);
                } else if (zdrop > 0) {
                    if (i - max_is[process_batch_id] >  mjs[process_batch_id] - max_js[process_batch_id]) {
                        //if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break;
                        if (maxs[process_batch_id] - ms[process_batch_id] - ((i - max_is[process_batch_id]) - ( mjs[process_batch_id] - max_js[process_batch_id])) * e_del <= zdrop) flag=1;
                    } else {
                        //if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break;
                        if (maxs[process_batch_id] - ms[process_batch_id]- (( mjs[process_batch_id] - max_js[process_batch_id]) - (i - max_is[process_batch_id])) * e_ins <= zdrop) flag=1;
                    }
                }
            }
            if(flag==0)
                break_flag=1;

            for(int process_batch_id = 0; process_batch_id<8; process_batch_id++)
            {
                int16_t j;
                for (j = begs[process_batch_id]; LIKELY(j < ends[process_batch_id]) && ehs[j][process_batch_id].h == 0 && ehs[j][process_batch_id].e == 0; ++j);
                begs[process_batch_id]=j;
                for (j = ends[process_batch_id]; LIKELY(j >= begs[process_batch_id]) && ehs[j][process_batch_id].h == 0 && ehs[j][process_batch_id].e == 0; --j);
                ends[process_batch_id] = j + 2 < qlens[process_batch_id]? j + 2 : qlens[process_batch_id];
          //  fprintf(stderr,"zdrop is %d\n",zdrop);
                // update beg and end for the next round
                //                    for (j = beg; LIKELY(j < end) && eh[j].h == 0 && eh[j].e == 0; ++j);
                //                    beg = j;
                //                    for (j = end; LIKELY(j >= beg) && eh[j].h == 0 && eh[j].e == 0; --j);
                //                    end = j + 2 < qlen? j + 2 : qlen;
                //beg = 0; end = qlen; // uncomment this line for debugging
            }
        }
        for(int process_batch_id = 0; process_batch_id<8; process_batch_id++)
        {
            g_qle[process_batch_id + grid_process_batch_idx*8] = max_js[process_batch_id]+1;
            g_tle[process_batch_id + grid_process_batch_idx*8] = max_is[process_batch_id]+1;
            g_gtle[process_batch_id + grid_process_batch_idx*8] = max_ies[process_batch_id]+1;
            g_gscore[process_batch_id + grid_process_batch_idx*8] = gscores[process_batch_id];
            g_max_off[process_batch_id + grid_process_batch_idx*8] = max_offs[process_batch_id];
            g_score[process_batch_id + grid_process_batch_idx*8] = maxs[process_batch_id];
            //db_hash_nxt_id++;
            
        }
        free(ehs);
    }
    free(qp_buff);
}
/**************/
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
            aligned_len = ((aligned_len+BATCHSIZE-1)/BATCHSIZE)*BATCHSIZE;
            uint64_t * tptr = &swlen_resized[i*BATCHSIZE];
            //last one in a batch
            global_batch_id = global_id_x*BATCHSIZE;
            for(int j=0; j<BATCHSIZE; j++)
            {
                hash_t* cur_hash_db = db_hash + i*BATCHSIZE+j;
                cur_hash_db->len = (uint32_t)tptr[j];
                cur_hash_db->alined = aligned_len;
                cur_hash_db->local_id_y = j;
                cur_hash_db->global_batch_id = global_batch_id;
            }
            global_id_x+=aligned_len;//point to next position
        }
        aligned_len =(uint32_t)swlen[size-1];//must be the last one
        aligned_len = ((aligned_len+BATCHSIZE-1)/BATCHSIZE)*BATCHSIZE;
        
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
    uint8_t* rdb_rev = malloc(global_id_x*BATCHSIZE);
    memset(rdb,0,sizeof(uint8_t)*global_id_x*BATCHSIZE);
    memset(rdb_rev,0,sizeof(uint8_t)*global_id_x*BATCHSIZE);
    //copy reference to db
    
    /***********************/
    //reference should be reverse
    //swlen_resize size is resize
    //construct RDB;
    for(int i=0; i<resize; i++)
    {
        swrst_t *sw = swrts+(swlen_resized[i]>>32);
        swseq_t *seq = sw->sw_seq;
        hash_t* rdb_hash_t = &db_hash[i];
        rdb_hash_t->rlen=seq->rlen;
        
        uint8_t* db_ptr = rdb+rdb_hash_t->global_batch_id+rdb_hash_t->local_id_y*rdb_hash_t->alined;
        
        memcpy(db_ptr,seq->ref,seq->rlen*sizeof(uint8_t));
    }
    for(int gride_batch_id = 0; gride_batch_id<resize_segs; gride_batch_id++)
    {
        hash_t* rdb_hash_t = &db_hash[gride_batch_id*BATCHSIZE];
        uint8_t* db_ptr = rdb + rdb_hash_t->global_batch_id;
        uint8_t* db_rev_ptr = rdb_rev + rdb_hash_t->global_batch_id;
        int x = BATCHSIZE;
        int y = rdb_hash_t->alined;
        transpose_8(db_ptr, db_rev_ptr, x, y);
    }
    //test
#ifdef DEBUG
    {
        hash_t* rdb_hash_t = &db_hash[1*BATCHSIZE];
        uint8_t* db_ptr = rdb + rdb_hash_t->global_batch_id;
        uint8_t* db_rev_ptr = rdb_rev + rdb_hash_t->global_batch_id;
        int size_x = BATCHSIZE;
        int size_y = rdb_hash_t->alined;
        for(int id_x=0; id_x<size_x; id_x++)
        {
            for(int id_y=0; id_y<size_y; id_y++)
            {
                assert(db_ptr[id_x*size_y+id_y]==db_rev_ptr[id_x+size_x*id_y]);
            }
        }
    }
#endif
    /***********************/
    //generate profile
    //qprofile should be aligned to 256bit (16x16)|(8X32)
    int16_t* qp_db = malloc(global_id_x*BATCHSIZE*g_m*sizeof(int16_t));//should be usingned ,change in the future
    memset(qp_db,(int16_t)-1,sizeof(int16_t)*global_id_x*BATCHSIZE*g_m);
    for(int i=0; i<resize; i++)
    {
        swrst_t *sw = swrts+(swlen_resized[i]>>32);
        swseq_t *seq = sw->sw_seq;
        hash_t* rdb_hash_t = &db_hash[i];
        
        rdb_hash_t->qlen=seq->qlen;
        int qlen =seq->qlen;
        int aligned = rdb_hash_t->alined;
        int16_t* qp = qp_db+g_m*(rdb_hash_t->global_batch_id+rdb_hash_t->local_id_y*rdb_hash_t->alined);
        const uint8_t* query =seq->query;
        for (int k = 0, m = 0; k < g_m; ++k) {
            const int8_t *p = g_mat[k];
            int j = 0;
            for (; j < qlen; ++j) qp[m++] = p[query[j]];
            for(;j<aligned; ++j) qp[m++]=p[5];
        }
    }
    /************************/
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
    for(int seg_idx=0; seg_idx<resize_segs;++seg_idx)
    {
        swlen_batch_id = swlen_resized+seg_idx* BATCHSIZE;
        db_hash_batch_id = db_hash+seg_idx*BATCHSIZE;
        int16_t g_qle[BATCHSIZE];
        int16_t g_tle[BATCHSIZE];
        int16_t g_gtle[BATCHSIZE];
        int16_t g_gscore[BATCHSIZE];
        int16_t g_max_off[BATCHSIZE];
        int16_t g_score[BATCHSIZE];
        int16_t g_h0[BATCHSIZE];// = sw->h0;
        swlen_nxt_id = swlen_batch_id;
        //threoretically wont work
        int batch_idx = 0;
        next_process = remain<BATCHSIZE?remain:BATCHSIZE;
#ifdef DEBUG
        fprintf(stderr,"remaining: %d\n",remain);
#endif
        remain-=BATCHSIZE;
        for(batch_idx=0; batch_idx<BATCHSIZE; batch_idx++)
        {
            g_h0[batch_idx] = 0;
        }
        for(batch_idx=0; batch_idx<next_process; batch_idx++)
        {
            swrst_t *sw = swrts+((*swlen_nxt_id)>>32);
            g_h0[batch_idx] = sw->h0;
            swlen_nxt_id++;
        }
        //process 16 query at a time
        batch_sw_core(db_hash_batch_id,
                            rdb_rev,
                            qp_db,
                            g_h0,//input
                            o_del,
                            e_del,
                            o_ins,
                            e_ins,
                            zdrop,
                           
                            g_qle,//result
                            g_tle,
                            g_gtle,
                            g_gscore,
                            g_max_off,
                            g_score
                           
                           );
        swlen_nxt_id = swlen_batch_id;
        for(batch_idx=0; batch_idx<next_process; batch_idx++)
        {
            swrst_t *sw = swrts+((*swlen_nxt_id)>>32);
            sw->qle = g_qle[batch_idx];
            sw->tle = g_tle[batch_idx];
            sw->gtle = g_gtle[batch_idx];
            sw->gscore = g_gscore[batch_idx];
            sw->max_off = g_max_off[batch_idx];
            sw->score = g_score[batch_idx];
            swlen_nxt_id++;
          //  if(swlen_nxt_id==swlen_end)break;
        }
    }
    
    
    free(qp_db);
    free(rdb);
    free(rdb_rev);
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
               //     if (m == 0) break;
                    if (m > max) {
                        max = m, max_i = i, max_j = mj;
                        max_off = max_off > abs(mj - i)? max_off : abs(mj - i);
                    } else if (zdrop > 0) {
//                        if (i - max_i > mj - max_j) {
//                            if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break;
//                        } else {
//                            if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break;
//                        }
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
    size_t process_sze = nread;
    fprintf(stderr,"now try %ld\n",process_sze);
    ksw_extend_batch2(nsrt, process_sze);
    
    //time
    fprintf(stderr, "[M::%s] Processed %ld reads in %.3f CPU sec, %.3f real sec\n", __func__, nread, cputime() - ctime, realtime() - rtime);

    swrst_t* rsrt;
    size_t rread = load(&rsrt,"sw_end_8000_0_2000.bin");
    assert(rread==nread);
    fprintf(stdout,"the latter result is correct\n");
    //printf("the result is %d\n", cmp(nsrt,rsrt,nread));
    uint8_t check = printdif(nsrt,rsrt,process_sze);
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
