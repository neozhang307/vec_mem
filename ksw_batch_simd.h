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
    int max_off;//64
    
    int qle,tle;//64
    int gtle, gscore;//64
    
    int h0;//64
    swseq_t* sw_seq;//64
}swrst_t;


