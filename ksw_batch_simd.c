#include"ksw_batch_simd.h"
#include <stdlib.h>

#include <assert.h>
#include <emmintrin.h>
#include<stdio.h>
#include<string.h>
void store(swrst_t* data, size_t size, const char* filename)
{
    FILE* output = fopen(filename,"wb+");
    fwrite(&size, sizeof(int), 1, output);
    
    for(int i=0; i<size; i++)
    {
        fwrite(data, sizeof(swrst_t), size, output);
    }
  
    for(int i=0; i<size; i++)
    {
        int qlen = data->sw_seq->qlen;
        const uint8_t* q = data->sw_seq->query;
        int rlen = data->sw_seq->qlen;
        const uint8_t* r = data->sw_seq->ref;
        fwrite(&qlen, sizeof(int), 1, output);
        fwrite(q, sizeof(uint8_t), qlen, output);
        fwrite(&rlen, sizeof(int), 1, output);
        fwrite(r, sizeof(uint8_t), rlen, output);
    }
    fclose(output);
}

size_t load(swrst_t** data, const char* filename)
{
    FILE* input = fopen(filename,"rb");
    if(input==NULL)
        return -1;
    int size;
    fread(&size, sizeof(int), 1, input);
    
    *data = malloc(sizeof(swrst_t)*size);
    swseq_t* seqs = malloc(sizeof(swseq_t)*size);
    for(int i=0; i<size; i++)
    {
        fread(*data, sizeof(swrst_t), size, input);
    }
    
    for(int i=0; i<size; i++)
    {
        int qlen;// = data->sw_seq->qlen;
        uint8_t* q ;// = data->sw_seq->query;
        int rlen;// = data->sw_seq->qlen;
        uint8_t* r ;//= data->sw_seq->ref;
        fread(&qlen, sizeof(int), 1, input);
        q = malloc(sizeof(uint8_t)* qlen);
        fread(q, sizeof(uint8_t), qlen, input);
        fread(&rlen, sizeof(int), 1, input);
        r = malloc(sizeof(uint8_t)* rlen);
        fread(r, sizeof(uint8_t), rlen, input);
        seqs[i].qlen=qlen;
        seqs[i].rlen=rlen;
        seqs[i].query = q;
        seqs[i].ref = r;
        data[0][i].sw_seq=seqs;
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
    free(data->sw_seq);
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
    fclose(confile);
}

void init(int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int zdrop)
{
    assert(m==5);
    memcpy(g_mat,mat,sizeof(int8_t)*25);
    g_o_del = o_del;
    g_e_del = e_del;
    g_o_ins = o_ins;
    g_e_del = e_del;
    g_zdrop = zdrop;
}