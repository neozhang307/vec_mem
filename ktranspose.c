//
//  main.cpp
//  MatrixTranspose
//
//  Created by zhanglingqi on 26/03/2017.
//  Copyright © 2017 zhanglingqi. All rights reserved.
//

#include <emmintrin.h>
//#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
//#define _MM_SHUFFLE(z, y, x, w) (z<<6) | (y<<4) | (x<<2) | w
//#define _128i_shuffle( m1, m2, z, y, x, w) _mm_shuffle_ps(m1,m2, _MM_SHUFFLE(z,y,x,w)
#define BATCHSIZE 8
#include<assert.h>
#define SHUFFLE(z, y, x, w) (w<<6) | (x<<4) | (y<<2) | z
inline void transpose_8x8_hwords (__m128i w0, __m128i w1,
                                  __m128i w2, __m128i w3,
                                  __m128i w4, __m128i w5,
                                  __m128i w6, __m128i w7,
                                  __m128i *r0, __m128i *r1,
                                  __m128i *r2, __m128i *r3,
                                  __m128i *r4, __m128i *r5,
                                  __m128i *r6, __m128i *r7)
{
    
    __m128i __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7;
    __m128i __tt0, __tt1, __tt2, __tt3, __tt4, __tt5, __tt6, __tt7;
    __t0 = _mm_unpacklo_epi16(w0, w1);
    __t1 = _mm_unpackhi_epi16(w0, w1);
    __t2 = _mm_unpacklo_epi16(w2, w3);
    __t3 = _mm_unpackhi_epi16(w2, w3);
    __t4 = _mm_unpacklo_epi16(w4, w5);
    __t5 = _mm_unpackhi_epi16(w4, w5);
    __t6 = _mm_unpacklo_epi16(w6, w7);
    __t7 = _mm_unpackhi_epi16(w6, w7);
    
    
    __tt0 = _mm_shuffle_ps(__t0,__t2,_MM_SHUFFLE(1,0,1,0));
    __tt1 = _mm_shuffle_ps(__t0,__t2,_MM_SHUFFLE(3,2,3,2));
    __tt2 = _mm_shuffle_ps(__t1,__t3,_MM_SHUFFLE(1,0,1,0));
    __tt3 = _mm_shuffle_ps(__t1,__t3,_MM_SHUFFLE(3,2,3,2));
    __tt4 = _mm_shuffle_ps(__t4,__t6,_MM_SHUFFLE(1,0,1,0));
    __tt5 = _mm_shuffle_ps(__t4,__t6,_MM_SHUFFLE(3,2,3,2));
    __tt6 = _mm_shuffle_ps(__t5,__t7,_MM_SHUFFLE(1,0,1,0));
    __tt7 = _mm_shuffle_ps(__t5,__t7,_MM_SHUFFLE(3,2,3,2));
    *r0 = _mm_shuffle_ps(__tt0, __tt4, _MM_SHUFFLE( 2, 0, 2, 0));
    *r1 = _mm_shuffle_ps(__tt0, __tt4, _MM_SHUFFLE( 3, 1, 3, 1));
    *r2 = _mm_shuffle_ps(__tt1, __tt5, _MM_SHUFFLE( 2, 0, 2, 0));
    *r3 = _mm_shuffle_ps(__tt1, __tt5, _MM_SHUFFLE( 3, 1, 3, 1));
    *r4 = _mm_shuffle_ps(__tt2, __tt6, _MM_SHUFFLE( 2, 0, 2, 0));
    *r5 = _mm_shuffle_ps(__tt2, __tt6, _MM_SHUFFLE( 3, 1, 3, 1));
    *r6 = _mm_shuffle_ps(__tt3, __tt7, _MM_SHUFFLE( 2, 0, 2, 0));
    *r7 = _mm_shuffle_ps(__tt3, __tt7, _MM_SHUFFLE( 3, 1, 3, 1));
}

inline void transpose_16x16_hwords (__m128i w0, __m128i w1,
                                    __m128i w2, __m128i w3,
                                    __m128i w4, __m128i w5,
                                    __m128i w6, __m128i w7,
                                    __m128i w8, __m128i w9,
                                    __m128i w10, __m128i w11,
                                    __m128i w12, __m128i w13,
                                    __m128i w14, __m128i w15,
                                    __m128i *r0, __m128i *r1,
                                    __m128i *r2, __m128i *r3,
                                    __m128i *r4, __m128i *r5,
                                    __m128i *r6, __m128i *r7,
                                    __m128i *r8, __m128i *r9,
                                    __m128i *r10, __m128i *r11,
                                    __m128i *r12, __m128i *r13,
                                    __m128i *r14, __m128i *r15)
{
    uint8_t* buff = (uint8_t*)malloc(sizeof(__m128i));
    
    __m128i __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7, __t8, __t9, __t10, __t11, __t12, __t13, __t14, __t15;
    __m128i __tt0, __tt1, __tt2, __tt3, __tt4, __tt5, __tt6, __tt7,__tt8, __tt9, __tt10, __tt11, __tt12, __tt13, __tt14, __tt15;
    //1
    __t0 = _mm_unpacklo_epi8(w0, w1);
    __t1 = _mm_unpackhi_epi8(w0, w1);
    __t2 = _mm_unpacklo_epi8(w2, w3);
    __t3 = _mm_unpackhi_epi8(w2, w3);
    __t4 = _mm_unpacklo_epi8(w4, w5);
    __t5 = _mm_unpackhi_epi8(w4, w5);
    __t6 = _mm_unpacklo_epi8(w6, w7);
    __t7 = _mm_unpackhi_epi8(w6, w7);
    __t8 = _mm_unpacklo_epi8(w8, w9);
    __t9 = _mm_unpackhi_epi8(w8, w9);
    __t10 = _mm_unpacklo_epi8(w10, w11);
    __t11 = _mm_unpackhi_epi8(w10, w11);
    __t12 = _mm_unpacklo_epi8(w12, w13);
    __t13 = _mm_unpackhi_epi8(w12, w13);
    __t14 = _mm_unpacklo_epi8(w14, w15);
    __t15 = _mm_unpackhi_epi8(w14, w15);
    _mm_store_si128((__m128i*)buff, __t0);
    //2
    __tt0 = _mm_unpacklo_epi16(__t0, __t2);
    __tt1 = _mm_unpackhi_epi16(__t0, __t2);
    __tt2 = _mm_unpacklo_epi16(__t1, __t3);
    __tt3 = _mm_unpackhi_epi16(__t1, __t3);
    __tt4 = _mm_unpacklo_epi16(__t4, __t6);
    __tt5 = _mm_unpackhi_epi16(__t4, __t6);
    __tt6 = _mm_unpacklo_epi16(__t5, __t7);
    __tt7 = _mm_unpackhi_epi16(__t5, __t7);
    __tt8 = _mm_unpacklo_epi16(__t8, __t10);
    __tt9 = _mm_unpackhi_epi16(__t8, __t10);
    __tt10 = _mm_unpacklo_epi16(__t9, __t11);
    __tt11 = _mm_unpackhi_epi16(__t9, __t11);
    __tt12 = _mm_unpacklo_epi16(__t12, __t14);
    __tt13 = _mm_unpackhi_epi16(__t12, __t14);
    __tt14 = _mm_unpacklo_epi16(__t13, __t15);
    __tt15 = _mm_unpackhi_epi16(__t13, __t15);
    //4
    __t0 = _mm_shuffle_ps(__tt0,__tt4,_MM_SHUFFLE(1,0,1,0));
    __t1 = _mm_shuffle_ps(__tt0,__tt4,_MM_SHUFFLE(3,2,3,2));
    __t2 = _mm_shuffle_ps(__tt1,__tt5,_MM_SHUFFLE(1,0,1,0));
    __t3 = _mm_shuffle_ps(__tt1,__tt5,_MM_SHUFFLE(3,2,3,2));
    __t4 = _mm_shuffle_ps(__tt2,__tt6,_MM_SHUFFLE(1,0,1,0));
    __t5 = _mm_shuffle_ps(__tt2,__tt6,_MM_SHUFFLE(3,2,3,2));
    __t6 = _mm_shuffle_ps(__tt3,__tt7,_MM_SHUFFLE(1,0,1,0));
    __t7 = _mm_shuffle_ps(__tt3,__tt7,_MM_SHUFFLE(3,2,3,2));
    __t8 = _mm_shuffle_ps(__tt8,__tt12,_MM_SHUFFLE(1,0,1,0));
    __t9 = _mm_shuffle_ps(__tt8,__tt12,_MM_SHUFFLE(3,2,3,2));
    __t10 = _mm_shuffle_ps(__tt9,__tt13,_MM_SHUFFLE(1,0,1,0));
    __t11 = _mm_shuffle_ps(__tt9,__tt13,_MM_SHUFFLE(3,2,3,2));
    __t12 = _mm_shuffle_ps(__tt10,__tt14,_MM_SHUFFLE(1,0,1,0));
    __t13 = _mm_shuffle_ps(__tt10,__tt14,_MM_SHUFFLE(3,2,3,2));
    __t14 = _mm_shuffle_ps(__tt11,__tt15,_MM_SHUFFLE(1,0,1,0));
    __t15 = _mm_shuffle_ps(__tt11,__tt15,_MM_SHUFFLE(3,2,3,2));
    
    //8
    *r0 = _mm_shuffle_ps(__t0,__t8,_MM_SHUFFLE(2,0,2,0));
    *r1 = _mm_shuffle_ps(__t0,__t8,_MM_SHUFFLE(3,1,3,1));
    *r2 = _mm_shuffle_ps(__t1,__t9,_MM_SHUFFLE(2,0,2,0));
    *r3 = _mm_shuffle_ps(__t1,__t9,_MM_SHUFFLE(3,1,3,1));
    *r4 = _mm_shuffle_ps(__t2,__t10,_MM_SHUFFLE(2,0,2,0));
    *r5 = _mm_shuffle_ps(__t2,__t10,_MM_SHUFFLE(3,1,3,1));
    *r6 = _mm_shuffle_ps(__t3,__t11,_MM_SHUFFLE(2,0,2,0));
    *r7 = _mm_shuffle_ps(__t3,__t11,_MM_SHUFFLE(3,1,3,1));
    *r8 = _mm_shuffle_ps(__t4,__t12,_MM_SHUFFLE(2,0,2,0));
    *r9 = _mm_shuffle_ps(__t4,__t12,_MM_SHUFFLE(3,1,3,1));
    *r10 = _mm_shuffle_ps(__t5,__t13,_MM_SHUFFLE(2,0,2,0));
    *r11 = _mm_shuffle_ps(__t5,__t13,_MM_SHUFFLE(3,1,3,1));
    *r12 = _mm_shuffle_ps(__t6,__t14,_MM_SHUFFLE(2,0,2,0));
    *r13 = _mm_shuffle_ps(__t6,__t14,_MM_SHUFFLE(3,1,3,1));
    *r14 = _mm_shuffle_ps(__t7,__t15,_MM_SHUFFLE(2,0,2,0));
    *r15 = _mm_shuffle_ps(__t7,__t15,_MM_SHUFFLE(3,1,3,1));
    
}

inline void transpose_4x4_dwords (__m128i w0, __m128i w1,
                                  __m128i w2, __m128i w3,
                                  __m128i *r0, __m128i *r1,
                                  __m128i *r2, __m128i *r3)
{
    
    __m128i x0 = _mm_shuffle_ps(w0, w1, _MM_SHUFFLE( 1, 0, 1,0)); // 0 1 4 5
    
    __m128i x1 = _mm_shuffle_ps (w0, w1,_MM_SHUFFLE(3, 2, 3, 2)); // 2 3 6 7
    
    __m128i x2 = _mm_shuffle_ps (w2, w3,_MM_SHUFFLE(1, 0, 1,0)); // 8 9 12 13
    
    __m128i x3 = _mm_shuffle_ps (w2, w3, _MM_SHUFFLE(3, 2, 3,2)); // 10 11 14 15
    
    *r0 = _mm_shuffle_ps (x0, x2,_MM_SHUFFLE( 2, 0, 2,0));//0 4 8 12
    *r1 = _mm_shuffle_ps (x0, x2,_MM_SHUFFLE( 3, 1, 3, 1));//1 5 9 13
    *r2 = _mm_shuffle_ps (x1, x3,_MM_SHUFFLE( 2, 0, 2,0));//2 6 10 14
    *r3 = _mm_shuffle_ps (x1, x3, _MM_SHUFFLE( 3, 1, 3, 1));//3 7 11 15
}

void transpose_16(int16_t*data_in, int16_t*data_out, int size_x, int size_y)
{
    __m128i*data_in_ori= (__m128i*)data_in;
    __m128i*data_out_ori=(__m128i*)data_out;
    __m128i*nxt_in = (__m128i*)data_in;
    __m128i*nxt_out = (__m128i*)data_out;
    
    int step_x = size_x/BATCHSIZE;
    assert(size_x%BATCHSIZE==0);
    
    int step_y = size_y/BATCHSIZE;
    assert(size_y%BATCHSIZE==0);
    
    for(int x=0; x<step_x; x++)
    {
        nxt_out = (__m128i*)data_out_ori + x;
        nxt_in = (__m128i*)data_in_ori + x*step_y*BATCHSIZE;
        for(int y=0; y<step_y; y++)
        {
            __m128i m0 = _mm_load_si128(nxt_in);
            
            __m128i m1 = _mm_load_si128(nxt_in+1*step_y);
            
            __m128i m2 = _mm_load_si128(nxt_in+2*step_y);
            
            __m128i m3 = _mm_load_si128(nxt_in+3*step_y);
            
            __m128i m4 = _mm_load_si128(nxt_in+4*step_y);
            
            __m128i m5 = _mm_load_si128(nxt_in+5*step_y);
            
            __m128i m6 = _mm_load_si128(nxt_in+6*step_y);
            __m128i m7 = _mm_load_si128(nxt_in+7*step_y);
            __m128i r0,r1,r2,r3,r4,r5,r6,r7;
            transpose_8x8_hwords(m0,m1,m2,m3,m4,m5,m6,m7,
                                 &r0,&r1,&r2,&r3, &r4,&r5,&r6,&r7);
            _mm_store_si128(nxt_out,r0);
            _mm_store_si128(nxt_out+1*step_x,r1);
            _mm_store_si128(nxt_out+2*step_x,r2);
            _mm_store_si128(nxt_out+3*step_x,r3);
            _mm_store_si128(nxt_out+4*step_x,r4);
            _mm_store_si128(nxt_out+5*step_x,r5);
            _mm_store_si128(nxt_out+6*step_x,r6);
            _mm_store_si128(nxt_out+7*step_x,r7);
            nxt_in +=1;
            nxt_out+=BATCHSIZE * step_x;
        }
    }
}
void transpose_8(uint8_t*data_in, uint8_t*data_out, int size_x, int size_y)
{
    __m128i*data_in_ori= (__m128i*)data_in;
    __m128i*data_out_ori=(__m128i*)data_out;
    __m128i*nxt_in = (__m128i*)data_in;
    __m128i*nxt_out = (__m128i*)data_out;
    
    int step_x = size_x/16;
    assert(size_x%16==0);
    
    int step_y = size_y/16;
    assert(size_y%16==0);
    
    for(int x=0; x<step_x; x++)
    {
        nxt_out = (__m128i*)data_out_ori + x;
        nxt_in = (__m128i*)data_in_ori + x*step_y*16;
        for(int y=0; y<step_y; y++)
        {
            __m128i m0 = _mm_load_si128(nxt_in);
            
            __m128i m1 = _mm_load_si128(nxt_in+1*step_y);
            
            __m128i m2 = _mm_load_si128(nxt_in+2*step_y);
            
            __m128i m3 = _mm_load_si128(nxt_in+3*step_y);
            
            __m128i m4 = _mm_load_si128(nxt_in+4*step_y);
            
            __m128i m5 = _mm_load_si128(nxt_in+5*step_y);
            
            __m128i m6 = _mm_load_si128(nxt_in+6*step_y);
            __m128i m7 = _mm_load_si128(nxt_in+7*step_y);
            __m128i m8 = _mm_load_si128(nxt_in+8*step_y);
            
            __m128i m9 = _mm_load_si128(nxt_in+9*step_y);
            
            __m128i m10 = _mm_load_si128(nxt_in+10*step_y);
            
            __m128i m11 = _mm_load_si128(nxt_in+11*step_y);
            
            __m128i m12 = _mm_load_si128(nxt_in+12*step_y);
            
            __m128i m13= _mm_load_si128(nxt_in+13*step_y);
            
            __m128i m14 = _mm_load_si128(nxt_in+14*step_y);
            __m128i m15= _mm_load_si128(nxt_in+15*step_y);
            __m128i r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15;
            transpose_16x16_hwords(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,
                                   &r0,&r1,&r2,&r3, &r4,&r5,&r6,&r7,&r8,&r9,&r10,&r11, &r12,&r13,&r14,&r15);
            _mm_store_si128(nxt_out,r0);
            _mm_store_si128(nxt_out+1*step_x,r1);
            _mm_store_si128(nxt_out+2*step_x,r2);
            _mm_store_si128(nxt_out+3*step_x,r3);
            _mm_store_si128(nxt_out+4*step_x,r4);
            _mm_store_si128(nxt_out+5*step_x,r5);
            _mm_store_si128(nxt_out+6*step_x,r6);
            _mm_store_si128(nxt_out+7*step_x,r7);
            _mm_store_si128(nxt_out+8*step_x,r8);
            _mm_store_si128(nxt_out+9*step_x,r9);
            _mm_store_si128(nxt_out+10*step_x,r10);
            _mm_store_si128(nxt_out+11*step_x,r11);
            _mm_store_si128(nxt_out+12*step_x,r12);
            _mm_store_si128(nxt_out+13*step_x,r13);
            _mm_store_si128(nxt_out+14*step_x,r14);
            _mm_store_si128(nxt_out+15*step_x,r15);
            nxt_in +=1;
            nxt_out+=16 * step_x;
        }
    }
}
//#define DEBUGTRANS
#ifdef DEBUGTRANS
int main(int argc, const char * argv[]) {
    // insert code here...
    int size_x = 16;
    int size_y = 16*5;
    uint8_t * space = (uint8_t*)malloc(size_x*size_y*sizeof(uint8_t));
    uint8_t * space_rst = (uint8_t*)malloc(size_x*size_y*sizeof(uint8_t));
    //uint16_t * space_buff = (uint16_t*)malloc(size_x*size_y);
    
    for(int i=0; i<size_y; i++)
    {
        for(int j=0; j<size_y; j++)
            space[i*size_y+j]=(i*size_y+j)%255;
    }
    for(int i=0; i<size_x; i++)
    {
        for(int j=0; j<size_y; j++)
            printf("%d ",space[i*size_y+j]);
        printf("\n");
    }
    printf("\n****\n");
    transpose(space,space_rst,size_x,size_y);
    
    for(int i=0; i<size_y; i++)
    {
        for(int j=0; j<size_x; j++)
            printf("%d ",space_rst[i*size_x+j]) ;
        printf("\n");
    }
    return 0;
}
#endif