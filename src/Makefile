NVCC = nvcc
CC = gcc

RS: matrix.o encode.o decode.o main.o
	$(NVCC) matrix.o encode.o decode.o main.o -o $@ 


.cu .o:
	$(NVCC) -o $@ -c $< 
#.c .o:
#	$(CC) -o $@ -c $<


#encode:
#	$(NVCC) -o encode.o encode.cu
clean:
	rm knn.o
test:
	$(CC) -o test-seq test-seq.c

include .depend

ChangeLog::
		hg log --style=changelog --only-branch=default "${top_srcdir}" > $@
