CC = /usr/bin/gcc
NVCC = nvcc
CCFLAGS = -g
# NVCCFLAGS = -G -arch=sm_35
NVCCFLAGS = -G -Xptxas -v

TARGET = RS-log-exp-V3 RS-log-exp-V0 RS-log-exp-V0-1 RS-log-exp-V1 RS-log-exp-V2
COMMON_OBJECTS = encode.o decode.o main.o
OBJECTS = encode.o decode.o main.o matrix.o matrix-V0.o matrix-V1.o matrix-V2.o
OBJDIR = ../bin
VPATH = $(OBJDIR)

OBJPROG = $(addprefix $(OBJDIR)/, $(TARGET))
all:$(OBJPROG)

TOTAL_OBJECTS := $(addprefix $(OBJDIR)/, $(OBJECTS))
$(TOTAL_OBJECTS): | $(OBJDIR)

$(OBJDIR)/%.o: %.cu	
	$(NVCC) -o $@ -c $< $(CCFLAGS) $(NVCCFLAGS)
$(OBJDIR)/%.o: %.c
	$(CC) -o $@ -c $< $(CCFLAGS)

$(OBJDIR)/RS-log-exp-V3: $(addprefix $(OBJDIR)/, $(COMMON_OBJECTS)) $(OBJDIR)/matrix.o
	$(NVCC) $^ -o $@ $(CCFLAGS) $(NVCCFLAGS)

$(OBJDIR)/RS-log-exp-V0: $(addprefix $(OBJDIR)/, $(COMMON_OBJECTS)) $(OBJDIR)/matrix-V0.o
	$(NVCC) $^ -o $@ $(CCFLAGS) $(NVCCFLAGS)

$(OBJDIR)/RS-log-exp-V0-1: $(addprefix $(OBJDIR)/, $(COMMON_OBJECTS)) $(OBJDIR)/matrix-V0-1.o
	$(NVCC) $^ -o $@ $(CCFLAGS) $(NVCCFLAGS)

$(OBJDIR)/RS-log-exp-V1: $(addprefix $(OBJDIR)/, $(COMMON_OBJECTS)) $(OBJDIR)/matrix-V1.o
	$(NVCC) $^ -o $@ $(CCFLAGS) $(NVCCFLAGS)

$(OBJDIR)/RS-log-exp-V2: $(addprefix $(OBJDIR)/, $(COMMON_OBJECTS)) $(OBJDIR)/matrix-V2.o
	$(NVCC) $^ -o $@ $(CCFLAGS) $(NVCCFLAGS)

$(OBJDIR):
	mkdir -p $(OBJDIR)

# test:
# 	$(CC) -o test-seq test-seq.c
CPU:
	$(CC) -o $(OBJDIR)/CPU-RS -lm -lrt $(CCFLAGS) cpu-rs.c

clean:
	rm $(OBJDIR)/*.o $(OBJPROG)
