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

void finalize_load(swrst_t*data)
{
    if(data==NULL)return;
    if(data->sw_seq!=NULL)
    {
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
        return -1;
    if(o_del!=g_o_del)
        return -2;
    if(e_del!=g_e_del)
        return -3;
    if(o_ins!=g_o_ins)
        return -5;
    if(e_ins!=g_e_ins)
        return -6;
    if(zdrop!=g_zdrop)
        return -7;
    for(int i=0; i<25; i++)
    {
        if(g_mat[0][i]!=mat[i])
            return -1*i*10;
    }
    return 1;
}
void printcmp(swrst_t* a, swrst_t*b, size_t size);
void printdif(swrst_t* a, swrst_t*b, size_t size)
{
    
   for(int i=0; i<size; i++)
    {
        printf("process %d of sw\n",i);
    
        swrst_t a_ = a[i];
        swrst_t b_ = b[i];
        if(a_.score!=b_.score)
        {
            printcmp(&a_,&b_,1);
            continue;
        }
        if(a_.qle!=b_.qle)
        {
            printcmp(&a_,&b_,1);
            continue;
        }if(a_.tle!=b_.tle)
        {
            printcmp(&a_,&b_,1);
            continue;
        }
        if(a_.gtle!=b_.gtle)
        {
            printcmp(&a_,&b_,1);
            continue;
        }
        if(a_.gscore!=b_.gscore)
        {
            printcmp(&a_,&b_,1);
            continue;
        }
        if(a_.max_off!=b_.max_off)
        {
            printcmp(&a_,&b_,1);
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
                break;
            }
        }
        assert(a_seq->rlen==b_seq->rlen);
        for(int i=0; i<a_seq->rlen; i++)
        {
            if(a_seq->ref[i]!=b_seq->ref[i])
            {
                printcmp(&a_,&b_,1);
                break;
            }
        }
    }
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
#ifdef SWBATCHDB
int main()
{
    load_config();
    swrst_t* nsrt;
    size_t nread = load(&nsrt,"sw_start_8000_0_2000.bin");
    
    for(int i=0; i<nread;++i)
    {
        swrst_t *sw = nsrt+i;
        swseq_t *seq = sw->sw_seq;
        if(seq->qlen!=0)
        {
            sw->score = ksw_extend2_mod(seq->qlen, seq->query, seq->rlen,seq->ref, 5, &g_mat[0][0], g_o_del, g_e_del, g_o_ins, g_e_ins, g_zdrop, sw->h0, &sw->qle, &sw->tle, &sw->gtle, &sw->gscore, &sw->max_off);
        }
    }
    swrst_t* rsrt;
    size_t rread = load(&rsrt,"sw_end_8000_0_2000.bin");
    assert(rread==nread);
    //printf("the result is %d\n", cmp(nsrt,rsrt,nread));
    printdif(nsrt,rsrt,nread);
    finalize_load(nsrt);
    finalize_load(rsrt);
    return 0;
}
#endif