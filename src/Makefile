NVCC = nvcc
CC = gcc
CCFLAGS = -g
NVCCFLAGS = -G

RS: matrix.o encode.o decode.o main.o
	$(NVCC) matrix.o encode.o decode.o main.o -o $@ 

%.o: %.cu
	$(NVCC) $(CCFLAGS) $(NVCCFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@

#.cu .o:
#	$(NVCC) -o $@ -c $< 

#.c .o:
#	$(CC) -o $@ -c $<


#encode:
#	$(NVCC) -o encode.o encode.cu
test:
	$(CC) -o test-seq test-seq.c
CPU:
	$(CC) -o CPU-RS -lm -lrt $(CCFLAGS) cpu-rs.c
clean:
	rm *.o && rm RS

include .depend

ChangeLog::
		hg log --style=changelog --only-branch=default "${top_srcdir}" > $@
