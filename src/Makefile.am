# AUTOMAKE_OPTIONS = subdir-objects

CCFLAGS = -g
NVCCFLAGS = -G

TARGET = RS
OBJECTS = matrix.o encode.o decode.o main.o

# OBJPROG = $(addprefix $(top_builddir)/, $(TARGET))

$(top_builddir)/RS: matrix.o encode.o decode.o main.o
	$(NVCC) $^ -o $@ 
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
# bin_PROGRAMS = RS
cuAutotools_SOURCES = matrix.cu encode.cu decode.cu main.c 

cuAutotools_CFLAGS  = $(CUDA_CFLAGS)
cuAutotools_LDADD  = $(CUDA_LIBS)
cuAutotools_LDFLAGS = 
