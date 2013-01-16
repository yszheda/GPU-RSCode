CCFLAGS = -g
NVCCFLAGS = -G

RS: matrix.o encode.o decode.o main.o
	$(NVCC) matrix.o encode.o decode.o main.o -o $@ 
.cu.o:
	$(NVCC) -o $@ -c $< $(CCFLAGS) $(NVCCFLAGS)
.c.o:
	$(CC) -o $@ -c $< $(CCFLAGS)
.cu.lo:
	$(top_srcdir)/cudalt.py $@ $(NVCC) $(NVCC_CFLAGS) --compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS)\" -c $<

#%.o: %.cu
#	$(NVCC) $(CCFLAGS) $(NVCCFLAGS) -c $< -o $@
#
#%.o: %.c
#	$(CC) $(CCFLAGS) -c $< -o $@

test:
	$(CC) -o test-seq test-seq.c
CPU:
	$(CC) -o CPU-RS -lm -lrt $(CCFLAGS) cpu-rs.c

#include .depend

CCLD=$(CC)
bin_PROGRAMS=RS
cuAutotools_SOURCES = matrix.cu encode.cu decode.cu main.c 

cuAutotools_CFLAGS  = $(CUDA_CFLAGS)
cuAutotools_LDADD  = $(CUDA_LIBS)
cuAutotools_LDFLAGS = 
