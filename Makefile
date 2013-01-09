NVCC = nvcc
CC = gcc
CCFLAGS = -g
NVCCFLAGS = -G
VPATH = src
#BIN_PATH = bin/

RS: matrix.o encode.o decode.o main.o
	$(NVCC) matrix.o encode.o decode.o main.o -o $@ 

#$(BIN_PATH)/%.o: %.cu
%.o: %.cu
	$(NVCC) $(CCFLAGS) $(NVCCFLAGS) -c $< -o $@

#$(BIN_PATH)/%.o: %.c
%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@

#.cu .o:
#	$(NVCC) -o $@ -c $< 

#.c .o:
#	$(CC) -o $@ -c $<

test:
	$(CC) -o test-seq $(VPATH)/test-seq.c
#	$(CC) -o $(BIN_PATH)/test-seq test-seq.c
CPU:
	$(CC) -o CPU-RS -lm -lrt $(CCFLAGS) $(VPATH)/cpu-rs.c
#	$(CC) -o $(BIN_PATH)/CPU-RS -lm -lrt $(CCFLAGS) cpu-rs.c
clean:
	rm *.o && rm RS

#include .depend
include $(VPATH)/.depend

ChangeLog::
		hg log --style=changelog --only-branch=default "${top_srcdir}" > $@
