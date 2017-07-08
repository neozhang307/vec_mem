gcc -O2 ssw.c -o ssw.o -fPIC -shared
gcc  -Wall -g -D  SWBATCHDB -O2 ssw.o ktranspose.o malloc_wrap.o ksw.o ksw_batch_simd.c -o ksw_batch.o
