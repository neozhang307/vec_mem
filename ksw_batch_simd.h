typedef struct
{
    int qlen;
    const uint8_t *query;
    int rlen;
    const uint8_t *ref;
}swseq_t;
typedef struct
{
    int score;
    int max_off;
    
    int qle,tle;
    int gtle, gscore;
    
    int h0;
    swseq_t* sw_seq;
}swrst_t;
