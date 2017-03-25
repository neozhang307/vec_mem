#include <stdint.h>
#include<stdio.h>
typedef struct
{
    uint32_t qlen;
    const uint8_t *query;
    uint32_t rlen;
    const uint8_t *ref;
}swseq_t;
typedef struct
{
    int score;
    int max_off;//64
    
    int qle,tle;//64
    int gtle, gscore;//64
    
    int h0;//64
    swseq_t* sw_seq;//64
}swrst_t;


void store(swrst_t* data, size_t size, const char* file);

size_t load(swrst_t** data, const char* file);

void finalize_load(swrst_t*data);

void load_config();
void store_config();
void init(int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int zdrop);

uint8_t printdif(swrst_t* a, swrst_t*b, size_t size);
int check_config(int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int zdrop);