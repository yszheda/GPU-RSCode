AUTOMAKE_OPTIONS = subdir-objects

CCFLAGS = -g
# NVCCFLAGS = -G
# NVCCFLAGS = -G -arch=sm_35
# NVCCFLAGS = -G -gencode=arch=compute_30,code=sm_30
NVCCFLAGS = -G -gencode=arch=compute_35,code=sm_35 --ptxas-options=-v
# NVCCFLAGS = -G --ptxas-options=-v

TARGET = RS RS-DFS-stream RS-BFS-stream
COMMON_OBJECTS = matrix.o decode.o main.o
OBJECTS = matrix.o encode.o encode-BFS-stream.o encode-DFS-stream.o decode.o main.o

# OBJPROG = $(addprefix $(top_builddir)/, $(TARGET))

$(top_builddir)/RS$(EXEEXT): matrix.o encode.o decode.o main.o
# $(top_builddir)/RS: matrix.o encode.o decode.o main.o
	$(NVCC) $^ -o $@ $(CCFLAGS) $(NVCCFLAGS) 

$(top_builddir)/RS-DFS-stream$(EXEEXT): $(COMMON_OBJECTS) encode-DFS-stream.o
# $(top_builddir)/RS-tiling-S1: $(COMMON_OBJECTS) matrix-S1.o
	$(NVCC) $^ -o $@ $(CCFLAGS) $(NVCCFLAGS) 

$(top_builddir)/RS-BFS-stream$(EXEEXT): $(COMMON_OBJECTS) encode-BFS-stream.o
# $(top_builddir)/RS-tiling-S1: $(COMMON_OBJECTS) matrix-S1.o
	$(NVCC) $^ -o $@ $(CCFLAGS) $(NVCCFLAGS) 

# $(OBJPROG): $(OBJECTS)
# 	$(NVCC) $^ -o $@ 
.cu.o:
	$(NVCC) -o $@ -c $< $(CCFLAGS) $(NVCCFLAGS)
.c.o:
	$(CC) -o $@ -c $< $(CCFLAGS)

# $(OBJDIR)/.lo: %.cu
# .cu.lo:
#	$(top_srcdir)/cudalt.py $@ $(NVCC) $(NVCC_CFLAGS) --compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS)\" -c $<


# test:
# 	$(CC) -o test-seq test-seq.c

CPU:
	$(CC) -o $(top_builddir)/CPU-RS -lm -lrt $(CCFLAGS) cpu-rs.c

# clean:
# 	rm *.o $(OBJPROG)
	
# include .depend

CCLD = $(CC)
bin_PROGRAMS = $(top_builddir)/RS $(top_builddir)/RS-DFS-stream $(top_builddir)/RS-BFS-stream

__top_builddir__RS_SOURCES = matrix.cu encode.cu decode.cu main.c 
__top_builddir__RS_BFS_stream_SOURCES = matrix.cu encode-BFS-stream.cu decode.cu main.c 
__top_builddir__RS_DFS_stream_SOURCES = matrix.cu encode-DFS-stream.cu decode.cu main.c 

# cuAutotools_CFLAGS  = $(CUDA_CFLAGS)
# cuAutotools_LDADD  = $(CUDA_LIBS)
# cuAutotools_LDFLAGS = 
