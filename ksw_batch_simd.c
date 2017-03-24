#include"ksw_batch_simd.h"
#include <stdlib.h>

#include <assert.h>
#include <emmintrin.h>
#include<stdio.h>
#include<string.h>
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