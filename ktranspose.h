#include <emmintrin.h>
//transpose for 16bit ellement, the size should be a multiply of 8
void transpose_i16(int16_t*data_in, int16_t*data_out, int size_x, int size_y);
//transpose for 16bit ellement, the size should be a multiply of 16
void transpose_u8(uint8_t*data_in, uint8_t*data_out, int size_x, int size_y);

void transpose_8x8_hwords (__m128i w0, __m128i w1,
                           __m128i w2, __m128i w3,
                           __m128i w4, __m128i w5,
                           __m128i w6, __m128i w7,
                           __m128i *r0, __m128i *r1,
                           __m128i *r2, __m128i *r3,
                           __m128i *r4, __m128i *r5,
                           __m128i *r6, __m128i *r7);

void transpose_16x16_hwords (__m128i w0, __m128i w1,
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
                             __m128i *r14, __m128i *r15);