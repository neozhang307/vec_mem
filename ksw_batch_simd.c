#include"ksw_batch_simd.h"
#include <stdlib.h>

#include <assert.h>
#include <emmintrin.h>
#include <tmmintrin.h>
#include<stdio.h>
#include<string.h>

#include"ksw.h"
#define ERR_NEXTLINE fprintf(stderr, "\n")
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
//static const int g_m = 5;
//static int8_t g_mat[5][5];
//static int g_o_del;
//static int g_e_del;
//static int g_o_ins;
//static int g_e_ins;
//static int g_zdrop;

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

void load_config( int8_t *mat, int *o_del, int *e_del, int *o_ins, int *e_ins, int *zdrop)
{
    FILE* confile = fopen("sw.config","rb");
    if(confile==NULL)
        return;
    fread(mat, sizeof(int8_t), 25, confile);
    fread(o_del, sizeof(int), 1, confile);
    fread(e_del, sizeof(int), 1, confile);
    fread(o_ins, sizeof(int), 1, confile);
    fread(e_ins, sizeof(int), 1, confile);
    fread(zdrop, sizeof(int), 1, confile);
    fclose(confile);
}

void store_config(const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int zdrop)
{
    FILE* confile = fopen("sw.config","wb+");
    if(confile==NULL)
        return;
    fwrite(mat, sizeof(int8_t), 25, confile);
    fwrite(&o_del, sizeof(int), 1, confile);
    fwrite(&e_del, sizeof(int), 1, confile);
    fwrite(&o_ins, sizeof(int), 1, confile);
    fwrite(&e_ins, sizeof(int), 1, confile);
    fwrite(&zdrop, sizeof(int), 1, confile);
    fclose(confile);
}

//void init(int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int zdrop)
//{
//    assert(m==5);
//    memcpy(g_mat,mat,sizeof(int8_t)*25);
////    g_o_del = o_del;
////    g_e_del = e_del;
////    g_o_ins = o_ins;
////    g_e_ins = e_ins;
////    g_zdrop = zdrop;
//}
//int check_config(int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int zdrop)
//{
//    if(m!=5)
//    {
//        return -1;
//    }
//    if(o_del!=g_o_del)
//    {
//        return -2;
//    }
//    if(e_del!=g_e_del)
//    {
//        return -3;
//    }
//    if(o_ins!=g_o_ins)
//    {
//        return -5;
//    }
//    if(e_ins!=g_e_ins)
//    {
//        return -6;
//    }
//    if(zdrop!=g_zdrop)
//    {
//        return -7;
//    }
//    for(int i=0; i<25; i++)
//    {
//        if(g_mat[0][i]!=mat[i])
//            return -1*i*10;
//    }
//    return 1;
//}
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

typedef struct{
    uint16_t len;
    uint16_t alined;
    uint16_t batch_max_len;
    uint8_t local_id_y;//idx inside batch
    size_t global_batch_id;//index of a batch in whole D
}packed_hash_t;
#define CHECK do{fprintf(stderr,"successfully process to line %d\n",__LINE__);}while(0)
#define CBUF(v_h0,h0s)    do{    _mm_store_si128((__m128i*)buffer,v_h0);\
        for(int process_batch_id=0; process_batch_id<8; process_batch_id++)\
        {\
            assert(h0s[process_batch_id]==buffer[process_batch_id]);\
        }\
}while(0)

void show(const char* str, int16_t* buffer, __m128i input)
{
    fprintf(stderr,"%s ", str);
    _mm_store_si128((__m128i*)buffer, input);
    for(int i=0; i<8; i++)
    {
        fprintf(stderr, " %d ",buffer[i]);
    }
    fprintf(stderr,"\n");
}

void batch_sw_core(packed_hash_t* ref_hash, packed_hash_t* que_hash,
                   uint8_t* rdb_rev,
                   int16_t* qp_db,
                   int16_t g_h0[BATCHSIZE],//input
                   
                   int m,
                   
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
    
  //  int16_t oe_del = o_del + e_del, oe_ins = o_ins + e_ins;
     __m128i v_zero, v_oe_del, v_e_del, v_oe_ins, v_e_ins,v_zdrop;
    v_zdrop = _mm_set1_epi16(zdrop);
    v_zero = _mm_set1_epi32(0);
    v_oe_del = _mm_set1_epi16(o_del + e_del);
    v_e_del = _mm_set1_epi16(e_del);
    v_oe_ins = _mm_set1_epi16(o_ins + e_ins);
    v_e_ins = _mm_set1_epi16(e_ins);

#define cmp_int16flag_change(ori_flag,v_zero,out_flag,tmp_h,tmp_l) do{\
(out_flag) = ori_flag;\
}while(0)
#define cmp_int16flag_change2(ori_flag,v_zero,out_flag,tmp_h,tmp_l) do{\
(tmp_h) = _mm_unpackhi_epi16(ori_flag, v_zero); \
(tmp_l) = _mm_unpacklo_epi16(v_zero, ori_flag);\
(tmp_h) = (__m128i)_mm_cmpneq_ps((__m128)tmp_h, (__m128)v_zero);\
(tmp_l) = (__m128i)_mm_cmpneq_ps((__m128)tmp_l, (__m128)v_zero);\
(out_flag) = _mm_packs_epi16(tmp_l, tmp_h);\
}while(0)
    //input would be modified
#define cmp_gen_result(cond,truecase,falsecase,tmp_out_true,tmp_out_false,out) do{\
(tmp_out_true)=(__m128i)_mm_and_ps((__m128)cond,(__m128)truecase);\
(tmp_out_false)=(__m128i)_mm_andnot_ps((__m128)cond,(__m128)falsecase);\
out = (__m128i)_mm_or_si128(tmp_out_true,tmp_out_false);\
}while(0)
    int que_align = que_hash->alined;
    
    size_t ref_batch_global_id = ref_hash->global_batch_id;
    size_t que_batch_global_id = que_hash->global_batch_id;
    
    int16_t* qp_buff = malloc(sizeof(int16_t)*PROCESSBATCH*que_align);
    memset(qp_buff,0,sizeof(int16_t)*PROCESSBATCH*que_align);
    int16_t* qp_buff_rev = malloc(sizeof(int16_t)*PROCESSBATCH*que_align);
    memset(qp_buff_rev,0,sizeof(int16_t)*PROCESSBATCH*que_align);
    
    uint16_t qlens[8];
    uint16_t maxqlen=0;
    uint16_t tlens[8];
    uint16_t maxtlen=0;
    __m128i v_qlen, v_tlen;//no smaller
    
  
   // int16_t begs[8], ends[8];

    
    //inreducable
    __m128i v_end;
    __m128i v_max, v_max_i, v_max_j, v_max_ie, v_gscore,v_max_off;
    //reducable
    __m128i v_h0 = _mm_set1_epi32(0);
    
    for(int grid_process_batch_idx=0; grid_process_batch_idx<BATCHSIZE/PROCESSBATCH;grid_process_batch_idx++)
    {
        const uint8_t *target_rev_batch =  rdb_rev+ref_batch_global_id;
        
        const int16_t *qp_batch = qp_db +m*que_batch_global_id;
        const int16_t *qp_batch_nxt = qp_batch + m*grid_process_batch_idx*8*que_align;
        
        //process 8 query at a time for int16_t
        v_h0=_mm_load_si128(((__m128i*)g_h0)+grid_process_batch_idx);
//        if(_mm_movemask_epi8(_mm_cmplt_ps((__m128)v_h0, (__m128)v_zero))!=0)//v_h0>0
//        {
//            
//        }
        assert(_mm_movemask_epi8(_mm_cmplt_epi16(v_h0, v_zero))==0);//all v_h0 > 0
        __m128i *v_hs = calloc(sizeof(__m128i)*(que_align+1),1);//can be smaller
        __m128i *v_es = calloc(sizeof(__m128i)*(que_align+1),1);//can be smaller
        
        v_hs[0]=v_h0;
        v_hs[1]=_mm_subs_epu16(v_h0, v_oe_ins);
        for(int j=2; j<=que_align; j++)
        {
            __m128i v_tmp = v_hs[j-1];
            
            v_tmp = _mm_subs_epu16(v_tmp, v_e_ins);
            
            if(!_mm_movemask_ps(_mm_cmpneq_ps((__m128)v_tmp, (__m128)v_zero)))
            {
                break;//when all equal to zero, break;
            }
            v_hs[j]=v_tmp;
        }
        /*********keeep************/
        for(int process_batch_id=0; process_batch_id<8; process_batch_id++)
        {
            // hash_t *tmphash =db_hash_nxt_id+process_batch_id+grid_process_batch_idx*8;
            packed_hash_t * tmp_quehash = que_hash+process_batch_id+grid_process_batch_idx*8;
            packed_hash_t * tmp_refhash = ref_hash+process_batch_id+grid_process_batch_idx*8;
            qlens[process_batch_id]=tmp_quehash->len;
            tlens[process_batch_id]=tmp_refhash->len;//tmphash->rlen;
        }
        maxqlen=que_hash->batch_max_len;
        maxtlen=ref_hash->batch_max_len;
        v_qlen=_mm_load_si128((__m128i*)qlens);
        
        v_tlen=_mm_load_si128((__m128i*)tlens);
        
        v_max = v_h0;
        v_max_i =_mm_set1_epi16(-1);
        v_max_j =_mm_set1_epi16(-1);
        v_max_ie =_mm_set1_epi16(-1);
        v_gscore =_mm_set1_epi16(-1);
        v_max_off = _mm_set1_epi32(0);
        /************************/
        /***********new*************/
        //reducable
        __m128i v_t, v_f, v_h1, v_m, v_mj;
        v_h1 = v_zero;
        //seems inreducable
        __m128i v_h_l, v_m_l, v_mj_l;
        int16_t min_beg, max_beg, min_end, max_end;
         /************************/

        v_end = v_qlen;

        __m128i tmplen = v_qlen;
        min_beg = 0;
        max_beg = 0;
        tmplen = v_end;
        __max_8(max_end, tmplen);
        tmplen = v_end;
        __max_8(min_end, tmplen);
        /************************/
        //MAIN SW
        for (int16_t i = 0; LIKELY(i < maxtlen) ; ++i) {
            __m128i v_tmp1;
            __m128i cond,cond2;
            __m128i truecase,falsecase;
            __m128i tmp_out_true,tmp_out_false;
            
            __m128i tmplen = v_end;
            __max_8(max_end, tmplen);
            //end possition should be updated
            /***********keep***********/
            uint8_t t_targets[8];
            memcpy(t_targets,target_rev_batch+i*BATCHSIZE + grid_process_batch_idx*8,8*sizeof(uint8_t));
            int16_t* qp_buff_nxt = qp_buff;
            for(int process_batch_id=0; process_batch_id<PROCESSBATCH; process_batch_id++)
            {
                
                int qp_ptr2 = process_batch_id*m*que_align+t_targets[process_batch_id] * que_align;
                
                memcpy(qp_buff_nxt, qp_batch_nxt+qp_ptr2, que_align*sizeof(int16_t));
                qp_buff_nxt+=que_align;
                
            }
            transpose_16(qp_buff,qp_buff_rev,PROCESSBATCH,que_align);
             __m128i v_i = _mm_set1_epi16(i);
            /***********************/
            uint16_t j;
            {

                //init
                v_f = v_zero;
                v_m = v_zero;
                v_mj = _mm_set1_epi16(-1);
                v_h_l = v_zero;
                v_m_l = v_zero;
                v_mj_l = v_zero;
                
                const  int16_t *q_rev = qp_buff_rev;
                // compute the first column
                if ( min_beg == 0) {
                    __m128i tval = _mm_set1_epi16((o_del + e_del * (i + 1)));
                    v_h1 = _mm_subs_epu16(v_h0,tval);
                } else v_h1=v_zero;
                
                //new reducable***********
                __m128i v_M;
                __m128i v_h;
                __m128i v_e;
                
                //processing a row
                    
                __m128i v_j = v_zero;
               
                //processing a row
               // for (j =  0; LIKELY(j < LOOP); ++j)
                for (j =  min_beg; LIKELY(j < max_end); ++j)
                {
                    v_j =_mm_set1_epi16(j);
                        // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
                        // Similar to SSE2-SW, cells are computed in the following order:
                        //   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
                        //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
                        //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
                    
                    v_M = v_hs[j];
                    v_e = v_es[j];

                    v_hs[j]=v_h1;
                    
                    cond =_mm_cmpeq_epi16(v_M,v_zero);
                    cmp_int16flag_change(cond, v_zero, cond, tmp_h, tmp_l);
                    __m128i tmp_qp = _mm_load_si128(((__m128i*)q_rev)+j);
                    falsecase = _mm_adds_epi16(v_M, tmp_qp);
                    cmp_gen_result(cond, v_zero, falsecase, tmp_out_true, tmp_out_false, v_M);
                   
                
                    v_h = _mm_max_epi16(v_M, v_e);
                    
                    v_h = _mm_max_epi16(v_h, v_f);

                    // save H(i,j) to h1 for the next column
                    v_h1 = v_h;
                   
                    cond = _mm_cmpgt_epi16(v_m, v_h);
                    truecase = v_mj;
                    falsecase = v_j;
                        
                    cmp_int16flag_change(cond,v_zero,cond,tmp_h,tmp_l);
                    cmp_gen_result(cond, v_mj, falsecase, tmp_out_true, tmp_out_false, v_mj);
                    cmp_gen_result(cond, v_m, v_h, tmp_out_true, tmp_out_false, v_m);

                    v_t=_mm_subs_epu16(v_M, v_oe_del);
                    
                    v_e = _mm_subs_epu16(v_e, v_e_del);
                    //condition+1
                    // computed E(i+1,j)
                    v_e = _mm_max_epi16(v_e, v_t);
                    v_es[j]=v_e;
                        
                        //saturation
                    v_t = _mm_subs_epu16(v_M, v_oe_ins);
                    v_f = _mm_subs_epu16(v_f, v_e_ins);
                    //condition+1
                    v_f = _mm_max_epi16(v_f, v_t);

                    //should think about it
                    cond =_mm_cmplt_epi16(v_j, v_end);
                    //redo unneccesary search
                    cmp_int16flag_change(cond, v_zero, cond, tmp_h, tmp_l);
                    cmp_gen_result(cond, v_h1, v_h_l, tmp_out_true, tmp_out_false, v_h_l);
                    cmp_gen_result(cond, v_m, v_m_l, tmp_out_true, tmp_out_false, v_m_l);
                    cmp_gen_result(cond, v_mj, v_mj_l, tmp_out_true, tmp_out_false, v_mj_l);
                    
                    
                }
                v_j =_mm_set1_epi16(j);

                v_m = v_m_l;
                v_mj = v_mj_l;
                v_h1 = v_h_l;
                
           //redo unneccesary search
                cond = _mm_cmplt_epi16(_mm_set1_epi16(i), v_tlen);
                cmp_int16flag_change(cond, v_zero, cond, tmp_h, tmp_l);
                cmp_gen_result(cond, v_h1, v_zero, tmp_out_true, tmp_out_false, v_h1);
                
                v_hs[j]=v_h1;
                v_es[j]=v_zero;
                
                v_j = _mm_min_epi16(v_j, v_end);
                cond = _mm_cmpeq_epi16(v_j, v_qlen);// when false no change   j==qlen?
                cmp_int16flag_change(cond, v_zero, cond, tmp_h, tmp_l);
                cond2 = _mm_cmpgt_epi16(v_gscore, v_h1);// when false no change//v_gscore<v_h1?
                // when true potentially change
                cmp_int16flag_change(cond2, v_zero, cond2, tmp_h, tmp_l);
                
                cond = (__m128i)_mm_andnot_ps((__m128)cond2, (__m128)cond);
    
                cmp_gen_result(cond, v_i, v_max_ie, tmp_out_true, tmp_out_false, v_max_ie);
                cmp_gen_result(cond, v_h1, v_gscore, tmp_out_true, tmp_out_false, v_gscore);
            }
                //if the search should terminated earlier?
            uint8_t flag = 0;
            
            //m==0 break
           
            if(!_mm_movemask_ps(_mm_cmpneq_ps((__m128)v_m, (__m128)v_zero)))
            {
                break;//break_flag=1;//when all equal to zero, break;
            }
            
            cond = _mm_cmpgt_epi16(v_m, v_max);// if (ms[process_batch_id] > maxs[process_batch_id])
            cmp_int16flag_change(cond, v_zero, cond, tmp_h, tmp_l);
            cmp_gen_result(cond, v_m, v_max, tmp_out_true, tmp_out_false, v_max);
            cmp_gen_result(cond, v_i, v_max_i, tmp_out_true, tmp_out_false, v_max_i);
            cmp_gen_result(cond, v_mj, v_max_j, tmp_out_true, tmp_out_false, v_max_j);
            //max_offs[process_batch_id] = max_offs[process_batch_id] > abs( mjs[process_batch_id] - i)? max_offs[process_batch_id] : abs( mjs[process_batch_id] - i);
            __m128i v_tmp_maxoff = _mm_abs_epi16( _mm_subs_epi16(v_mj, v_i));
            v_tmp_maxoff = _mm_max_epi16(v_tmp_maxoff, v_max_off);
            cmp_gen_result(cond, v_tmp_maxoff, v_max_off, tmp_out_true, tmp_out_false, v_max_off);

            flag=1;

            if(zdrop>0&&!_mm_movemask_epi8(cond))//all false
            {
                __m128i v_tmp2,v_tmp3,v_tmp4,v_tmp5,v_tmp6,v_tmp7,v_tmp8,cond3;
                flag=0;
                v_tmp1 = _mm_subs_epi16(v_i, v_max_i);//i - max_is[process_batch_id]
                v_tmp2 = _mm_subs_epi16(v_mj, v_max_j);//mjs[process_batch_id] - max_js[process_batch_id]
                cond = _mm_cmpgt_epi16(v_tmp1, v_tmp2);
                cmp_int16flag_change(cond, v_zero, cond, tmp_h, tmp_l);
                v_tmp3 = _mm_subs_epi16(v_tmp1, v_tmp2);//(i - max_is[process_batch_id]) - (mjs[process_batch_id] - max_js[process_batch_id])
                v_tmp4 = _mm_mulhrs_epi16(v_tmp3, v_e_del);//(i - max_is[process_batch_id]) - (mjs[process_batch_id] - max_js[process_batch_id])*e_del
                v_tmp5 = _mm_mulhrs_epi16(v_tmp3, v_e_ins);//(i - max_is[process_batch_id]) - (mjs[process_batch_id] - max_js[process_batch_id])*e_ins
                v_tmp6 = _mm_subs_epi16(v_max, v_m);//maxs[process_batch_id] - ms[process_batch_id]
                
                v_tmp7 = _mm_subs_epi16(v_tmp6, v_tmp4);
                v_tmp8 = _mm_adds_epi16(v_tmp6, v_tmp5);
                cond2 = _mm_cmplt_epi16(v_tmp7, v_zdrop);
                cmp_int16flag_change(cond2, v_zero, cond2, tmp_h, tmp_l);
         
                cond3 = _mm_cmplt_epi16(v_tmp8, v_zdrop);
                cmp_int16flag_change(cond3, v_zero, cond3, tmp_h, tmp_l);
               
                cond2 = _mm_and_si128(cond, cond2);
                cond3 = _mm_andnot_si128(cond, cond3);
                
                cond = _mm_or_si128(cond2, cond3);
                if(_mm_movemask_epi8(cond))flag=1;
                
            }

            v_tmp1 = v_end;
            __min_8(min_end,v_tmp1);
            v_tmp1 = v_end;
            __max_8(max_end,v_tmp1);
            
            
            for(j=min_beg; LIKELY(j<min_end); j++)
            {
                if(_mm_movemask_ps(_mm_cmpneq_ps((__m128)v_hs[j],(__m128)v_zero)))break;//any one not zero break
                if(_mm_movemask_ps(_mm_cmpneq_ps((__m128)v_es[j],(__m128)v_zero)))break;//any one not zero break
            }
             min_beg = j;
            for(; LIKELY(j<min_end); j++)
            {
                if(!_mm_movemask_ps(_mm_cmpeq_ps((__m128)v_hs[j],(__m128)v_zero)))break;//all not zero break
                if(!_mm_movemask_ps(_mm_cmpeq_ps((__m128)v_es[j],(__m128)v_zero)))break;//all not zero break
            }
             max_beg = j;
            
            for(j=max_end; LIKELY(j>max_beg); j--)
            {
                if(_mm_movemask_ps(_mm_cmpneq_ps((__m128)v_hs[j],(__m128)v_zero)))break;//any one not zero break
                if(_mm_movemask_ps(_mm_cmpneq_ps((__m128)v_es[j],(__m128)v_zero)))break;//any one not zero break
            }
            max_end = j;
            for(;LIKELY(j>max_beg); j--)
            {
                if(!_mm_movemask_ps(_mm_cmpeq_ps((__m128)v_hs[j],(__m128)v_zero)))break;//all not zero break
                if(!_mm_movemask_ps(_mm_cmpeq_ps((__m128)v_es[j],(__m128)v_zero)))break;//all not zero break
            }
            min_end = j;
            v_tmp1 = _mm_set1_epi16(max_end+2);
            v_end = _mm_min_epi16(v_tmp1, v_qlen);
        }
        __m128i v_tmp1;
        __m128i v_one = _mm_set1_epi16(1);
        v_tmp1 = _mm_add_epi16(v_max_j, v_one);
        _mm_store_si128(((__m128i*)g_qle)+grid_process_batch_idx,v_tmp1);
        v_tmp1 = _mm_add_epi16(v_max_i, v_one);
        _mm_store_si128(((__m128i*)g_tle)+grid_process_batch_idx,v_tmp1);
        v_tmp1 = _mm_add_epi16(v_max_ie, v_one);
        _mm_store_si128(((__m128i*)g_gtle)+grid_process_batch_idx,v_tmp1);
       
        _mm_store_si128(((__m128i*)g_gscore)+grid_process_batch_idx,v_gscore);
        _mm_store_si128(((__m128i*)g_max_off)+grid_process_batch_idx,v_max_off);
        _mm_store_si128(((__m128i*)g_score)+grid_process_batch_idx,v_max);
        
        free(v_hs);
        free(v_es);

    }
    free(qp_buff);
    free(qp_buff_rev);

}


void batch_sw_core2(packed_hash_t* ref_hash, packed_hash_t* que_hash,
                   uint8_t* rdb_rev,
                   int16_t* qp_db,
                   int16_t g_h0[BATCHSIZE],//input
                   
                    int m,
                    
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
    
    int16_t oe_del = o_del + e_del, oe_ins = o_ins + e_ins;
    int que_align = que_hash->alined;

    size_t ref_batch_global_id = ref_hash->global_batch_id;
    size_t que_batch_global_id = que_hash->global_batch_id;
    
    int16_t* qp_buff = malloc(sizeof(int16_t)*PROCESSBATCH*que_align);
    memset(qp_buff,0,sizeof(int16_t)*PROCESSBATCH*que_align);
    int16_t* qp_buff_rev = malloc(sizeof(int16_t)*PROCESSBATCH*que_align);
    memset(qp_buff_rev,0,sizeof(int16_t)*PROCESSBATCH*que_align);
    
    for(int grid_process_batch_idx=0; grid_process_batch_idx<BATCHSIZE/PROCESSBATCH;grid_process_batch_idx++)
    {
        const uint8_t *target_rev_batch =  rdb_rev+ref_batch_global_id;
 
        const int16_t *qp_batch = qp_db +m*que_batch_global_id;
        const int16_t *qp_batch_nxt = qp_batch + m*grid_process_batch_idx*8*que_align;
        
        //process 8 query at a time for int16_t
        uint16_t qlens[8];
        uint16_t maxqlen=0;
        uint16_t tlens[8];
        uint16_t maxtlen=0;
        
        eh_m **ehs=malloc(sizeof(eh_m*)*(que_align + 1));
        for(int i=0; i<que_align+1;i++)
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
            for (int j = 2; j <= que_align && ehs[j-1][process_batch_id].h > e_ins; ++j)
                ehs[j][process_batch_id].h = ehs[j-1][process_batch_id].h - e_ins;
        }
        for(int process_batch_id=0; process_batch_id<8; process_batch_id++)
        {
           // hash_t *tmphash =db_hash_nxt_id+process_batch_id+grid_process_batch_idx*8;
            packed_hash_t * tmp_quehash = que_hash+process_batch_id+grid_process_batch_idx*8;
            packed_hash_t * tmp_refhash = ref_hash+process_batch_id+grid_process_batch_idx*8;
            qlens[process_batch_id]=tmp_quehash->len;
            maxqlen=maxqlen>tmp_quehash->len?maxqlen:tmp_quehash->len;
            
            tlens[process_batch_id]=tmp_refhash->len;//tmphash->rlen;
            maxtlen=tmp_refhash->batch_max_len;//maxtlen>tmphash->rlen?maxtlen:tmphash->rlen;
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
            
 
            uint8_t t_targets[8];
            memcpy(t_targets,target_rev_batch+i*BATCHSIZE + grid_process_batch_idx*8,8*sizeof(uint8_t));
            int16_t* qp_buff_nxt = qp_buff;
            for(int process_batch_id=0; process_batch_id<PROCESSBATCH; process_batch_id++)
            {
             
                int qp_ptr2 = process_batch_id*m*que_align+t_targets[process_batch_id] * que_align;
                
                memcpy(qp_buff_nxt, qp_batch_nxt+qp_ptr2, que_align*sizeof(int16_t));
                qp_buff_nxt+=que_align;

            }
            transpose_16(qp_buff,qp_buff_rev,PROCESSBATCH,que_align);
            //transpose_16(qp_buff,qp_buff_rev,PROCESSBATCH,align_end);
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
                

                /***********************/

                const  int16_t *q_rev = qp_buff_rev;
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
                    
                    M = M? M + q_rev[process_batch_id+j*PROCESSBATCH] : 0;// separating H and M to disallow a cigar like "100M3I3D20M"
                    
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
        
        for(int i=0; i<que_align+1;i++)
        {
            free(ehs[i]);//=calloc(sizeof(eh_m)*8,1);
        }
        free(ehs);
    }
    free(qp_buff);
    free(qp_buff_rev);
}
/**************/
#include<time.h>
void ksw_extend_batch2(swrst_t* swrts, uint32_t size, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int zdrop)
{
#ifdef SWBATCHDB
    clock_t m_start,m_end;

    m_start = clock();
#endif
    //sort
    assert(m==5);
    assert(size>=0);
    uint64_t* swlen = malloc(sizeof(int64_t)*size);//should record qlen rlen
    int j=0;
    for(uint32_t i=0; i<size; ++i)
    {
        uint64_t tval=(uint64_t)i<<32;
        swseq_t* tseq = swrts[i].sw_seq;
        uint32_t qlen = tseq->qlen;
        if(qlen==0)continue;
 //       uint32_t rlen = tseq->rlen;
//        uint32_t mx = qlen>rlen?qlen:rlen;
        uint32_t mx=qlen;
        tval|=mx;
        swlen[j++]=tval;
    }
//    ks_introsort(uint64_t,size,swlen);
#ifdef DEBUG_SW
    assert((uint32_t)swlen[0]==0);
#endif
    int ptr = 0;
    
//    int threashold = 0;
//    while((uint32_t)swlen[ptr]==0&&ptr<size) ptr++;
//    
//    int none_zero = ptr;
//    while((uint32_t)swlen[ptr]<threashold&&ptr<size) ptr++;
//#ifdef SWBATCHDB
//    fprintf(stderr,"jump from %d to %d, while total size is%d\n",none_zero,ptr,size);
//#endif
//    //    //sinple use sinple way
//    for(int i=none_zero; i<ptr; i++)
//    {
//        int idx = swlen[i]>>32;
//        swrst_t *sw = swrts+idx;
//        swseq_t *seq = sw->sw_seq;
//        sw->score = ksw_extend2_mod(seq->qlen, seq->query, seq->rlen,seq->ref, m, mat, o_del, e_del, o_ins, e_ins, zdrop, sw->h0, &sw->qle, &sw->tle, &sw->gtle, &sw->gscore, &sw->max_off);
//    }
//
    size = j;
    //recompute space size
    //uint8_t batch = BATCHSIZE;//16 | 8
    
    uint32_t resize = size-ptr;
    uint32_t resize_segs = (resize+BATCHSIZE-1)/BATCHSIZE;//skip zero ones
    uint32_t aligned_resize = resize_segs*BATCHSIZE;
    uint64_t *swlen_resized = swlen+ptr;
    
    packed_hash_t* ref_hash = malloc(sizeof(packed_hash_t)*aligned_resize);
    memset(ref_hash,0,sizeof(packed_hash_t)*aligned_resize);
    
    packed_hash_t* que_hash = malloc(sizeof(packed_hash_t)*aligned_resize);
    memset(que_hash,0,sizeof(packed_hash_t)*aligned_resize);

    
    //init
    int ref_global_id_x=0;
    {
        int ref_aligned_len=0;
        int ref_global_batch_id = 0;
        int rlen_max = 0;
        for(int i=0; i<resize_segs; i++)
        {
            for(int j=0,k=i*BATCHSIZE; k<resize&&j<BATCHSIZE; j++,k++)
            {
                swrst_t *sw = swrts+(swlen_resized[k]>>32);
                swseq_t *seq = sw->sw_seq;
                int rlen = seq->rlen;
                rlen_max = rlen_max>rlen?rlen_max:rlen;
            }
            ref_aligned_len = ((rlen_max+BATCHSIZE-1)/BATCHSIZE)*BATCHSIZE;
            ref_global_batch_id = ref_global_id_x*BATCHSIZE;
            for(int j=0; j<BATCHSIZE; j++)
            {
                packed_hash_t *cur_ref_hash = ref_hash+i*BATCHSIZE+j;
                cur_ref_hash->local_id_y=j;
                cur_ref_hash->global_batch_id=ref_global_batch_id;
                cur_ref_hash->alined=ref_aligned_len;
                cur_ref_hash->batch_max_len=rlen_max;
            }
            ref_global_id_x+=ref_aligned_len;
        }
    }
    
    int que_global_id_x=0;
    {
        int que_aligned_len=0;
        int que_global_batch_id = 0;
        int qlen_max = 0;
        for(int i=0; i<resize_segs; i++)
        {
            for(int j=0,k=i*BATCHSIZE; k<resize&&j<BATCHSIZE; j++,k++)
            {
                swrst_t *sw = swrts+(swlen_resized[k]>>32);
                swseq_t *seq = sw->sw_seq;
                int qlen = seq->qlen;
                qlen_max = qlen_max>qlen?qlen_max:qlen;
            }
            que_aligned_len = ((qlen_max+BATCHSIZE-1)/BATCHSIZE)*BATCHSIZE;
            que_global_batch_id = que_global_id_x*BATCHSIZE;
            
            
            for(int j=0; j<BATCHSIZE; j++)
            {
                packed_hash_t *cur_que_hash = que_hash+i*BATCHSIZE+j;
                cur_que_hash->local_id_y=j;
                cur_que_hash->global_batch_id=que_global_batch_id;
                cur_que_hash->alined=que_aligned_len;
                cur_que_hash->batch_max_len=qlen_max;
//                fprintf(stderr,"global_batch_id %ld, local_id_y %d, align %d , batch_max_len %d \n", cur_que_hash->global_batch_id,cur_que_hash->local_id_y,cur_que_hash->alined, cur_que_hash->batch_max_len);
            }
            que_global_id_x+=que_aligned_len;
        }
    }
    
    

    uint8_t* rdb = malloc(ref_global_id_x*BATCHSIZE);
    uint8_t* rdb_rev = malloc(ref_global_id_x*BATCHSIZE);
    memset(rdb,0,sizeof(uint8_t)*ref_global_id_x*BATCHSIZE);
    memset(rdb_rev,0,sizeof(uint8_t)*ref_global_id_x*BATCHSIZE);
    
    //copy reference to db
    
    /***********************/
    //reference should be reverse
    //swlen_resize size is resize
    //construct RDB;
#ifdef SWBATCHDB
    clock_t start,end;
    start = clock();
#endif
    for(int i=0; i<resize; i++)
    {
        swrst_t *sw = swrts+(swlen_resized[i]>>32);
        swseq_t *seq = sw->sw_seq;
        packed_hash_t* ref_hash_t = &ref_hash[i];
        ref_hash_t->len=seq->rlen;
        
        uint8_t* db_ptr = rdb+ref_hash_t->global_batch_id+ref_hash_t->local_id_y*ref_hash_t->alined;
        memcpy(db_ptr,seq->ref,seq->rlen*sizeof(uint8_t));
    }
#ifdef SWBATCHDB
    end = clock();
    fprintf(stderr,"ref copy time:%f\n",(float)(end - start) / CLOCKS_PER_SEC);
#endif
#ifdef DEBUG_SW
    for(int i=0; i<resize; i++)
    {
        swrst_t *sw = swrts+(swlen_resized[i]>>32);
        swseq_t *seq = sw->sw_seq;
        packed_hash_t* ref_hash_t = &ref_hash[i];
        uint8_t* db_ptr = rdb+ref_hash_t->global_batch_id+ref_hash_t->local_id_y*ref_hash_t->alined;
        for(int j=0; j<seq->rlen; j++)
        {
        //    fprintf(stderr,"%d:%d",seq->ref[j],db_ptr2[j]);
            assert(seq->ref[j]==db_ptr[j]);
        }
    }
#endif
#ifdef SWBATCHDB
    start = clock();
#endif
    for(int gride_batch_id = 0; gride_batch_id<resize_segs; gride_batch_id++)
    {
        packed_hash_t* ref_hash_t = &ref_hash[gride_batch_id*BATCHSIZE];
        uint8_t* db_ptr = rdb + ref_hash_t->global_batch_id;
        uint8_t* db_rev_ptr = rdb_rev + ref_hash_t->global_batch_id;
        int x = BATCHSIZE;
        int y = ref_hash_t->alined;
        transpose_8(db_ptr, db_rev_ptr, x, y);
    }
#ifdef SWBATCHDB
    end = clock();
    fprintf(stderr,"ref transpose time:%f\n",(float)(end - start) / CLOCKS_PER_SEC);
#endif
 #ifdef DEBUG_SW
    {
        packed_hash_t* rdb_hash_t = &ref_hash[1*BATCHSIZE];
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
    int16_t* qp_db = malloc(que_global_id_x*BATCHSIZE*m*sizeof(int16_t));//should be usingned ,change in the future
    memset(qp_db,(int16_t)-1,sizeof(int16_t)*que_global_id_x*BATCHSIZE*m);
#ifdef SWBATCHDB
    start = clock();
#endif
    for(int i=0; i<resize; i++)
    {

        swrst_t *sw = swrts+(swlen_resized[i]>>32);
        swseq_t *seq = sw->sw_seq;
        packed_hash_t* que_hash_cur_ptr = &que_hash[i];

        que_hash_cur_ptr->len=seq->qlen;
        //fprintf(stderr, "%d:%d ",seq->qlen,que_hash_cur_ptr->len);
        
        int qlen =seq->qlen;
        //assert(qlen==1);
        int aligned = que_hash_cur_ptr->alined;
        
        //fprintf(stderr,"global_batch_id %ld, local_id_y %d, align %d \n",que_hash_t->global_batch_id,que_hash_t->local_id_y,que_hash_t->alined);
        int16_t* qp2 = qp_db+m*(que_hash_cur_ptr->global_batch_id+que_hash_cur_ptr->local_id_y*que_hash_cur_ptr->alined);
        
        //fprintf(stderr,"check%d\n",qp2[0]);
        const uint8_t* query =seq->query;
        for (int k = 0, l = 0; k < m; ++k) {
            const int8_t *p = &mat[k*m];
            int j = 0;
            for (; j < qlen; ++j) qp2[l++] = p[query[j]];
            for(;j<aligned; ++j) qp2[l++]=p[5];
        }
    }
#ifdef SWBATCHDB
    end = clock();
    fprintf(stderr,"qp execution time:%f\n",(float)(end - start) / CLOCKS_PER_SEC);
#endif
//    
    
    /************************/
  //  for(int i=ptr; i<size;++i)
    

 
    
    uint64_t* swlen_batch_id = swlen_resized;//resize
    packed_hash_t * ref_hash_batch_ptr = ref_hash;
    packed_hash_t * que_hash_batch_ptr = que_hash;
    uint64_t* swlen_nxt_id = swlen_resized;//resize
    uint32_t remain = resize;
    uint32_t next_process = BATCHSIZE;
    //main
#ifdef SWBATCHDB
    start = clock();
#endif
    for(int seg_idx=0; seg_idx<resize_segs;++seg_idx)
    {
        swlen_batch_id = swlen_resized+seg_idx* BATCHSIZE;
        ref_hash_batch_ptr = ref_hash + seg_idx*BATCHSIZE;
        que_hash_batch_ptr = que_hash + seg_idx*BATCHSIZE;
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
        batch_sw_core(ref_hash_batch_ptr,que_hash_batch_ptr,
                            rdb_rev,
                            qp_db,
                            g_h0,//input
                      
                            m,
                      
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
        }
    }
#ifdef SWBATCHDB
    end = clock();
    fprintf(stderr,"main execution time:%f\n",(float)(end - start) / CLOCKS_PER_SEC);
#endif
    free(ref_hash);
    free(que_hash);
    free(rdb);
    free(rdb_rev);
    free(swlen);
    free(qp_db);
#ifdef SWBATCHDB
    m_end = clock();
    fprintf(stderr,"total execution time:%f\n",(float)(m_end - m_start) / CLOCKS_PER_SEC);
#endif
}



void ksw_extend_batch(swrst_t* swrts, size_t size, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int zdrop)
{
    //sort
    assert(m==5);
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
         
            int h0 = sw->h0;

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
//#define SWBATCHDB
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
     const int g_m = 5;
     int8_t g_mat[5][5];
     int g_o_del;
     int g_e_del;
     int g_o_ins;
     int g_e_ins;
     int g_zdrop;
    
    load_config(&g_mat[0][0], &g_o_del, &g_e_del, &g_o_ins, &g_e_ins, &g_zdrop);
    swrst_t* nsrt;
    size_t nread = load(&nsrt,"sw_start_8000_0_2000.bin");
    
    //time
    double ctime = cputime();
    double rtime = realtime();
    size_t process_sze = nread;
    fprintf(stderr,"now try %ld\n",process_sze);
    ksw_extend_batch2(nsrt, process_sze, g_m, g_mat[0], g_o_del, g_e_del, g_o_ins, g_e_ins,g_zdrop);
    
    //time
    fprintf(stderr, "[M::%s] Processed %ld reads in %.3f CPU sec, %.3f real sec\n", __func__, nread, cputime() - ctime, realtime() - rtime);

//    swrst_t* rsrt;
//    size_t rread = load(&rsrt,"sw_end_8000_0_2000.bin");
//    assert(rread==nread);
//    fprintf(stdout,"the latter result is correct\n");
// //   printf("the result is %d\n", cmp(nsrt,rsrt,nread));
//    uint8_t check = printdif(nsrt,rsrt,process_sze);
//   if(check==0)
//    {
//        fprintf(stderr,"the check result is correct\n");
//        fprintf(stderr,"which is 0x%02x\n",check);
//    }
//    else
//    {
//        fprintf(stderr,"the check result is incorrect\n");
//        fprintf(stderr,"which is 0x%02x\n",check);
//    }
    finalize_load(nsrt,nread);
//    finalize_load(rsrt,rread);
    return 0;
}
#endif
