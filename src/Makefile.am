AUTOMAKE_OPTIONS = subdir-objects

CCFLAGS = -g
NVCCFLAGS = -G --ptxas-options=-v

TARGET = RS RS-tiling-S1 RS-tiling-S2 RS-tiling-S3
COMMON_OBJECTS = encode.o decode.o main.o
OBJECTS = matrix.o encode.o decode.o main.o

# OBJPROG = $(addprefix $(top_builddir)/, $(TARGET))

$(top_builddir)/RS$(EXEEXT): matrix.o encode.o decode.o main.o
# $(top_builddir)/RS: matrix.o encode.o decode.o main.o
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
bin_PROGRAMS = $(top_builddir)/RS

__top_builddir__RS_SOURCES = matrix.cu encode.cu decode.cu main.c 

# cuAutotools_CFLAGS  = $(CUDA_CFLAGS)
# cuAutotools_LDADD  = $(CUDA_LIBS)
# cuAutotools_LDFLAGS = 
