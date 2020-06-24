PLASMA_DIR= ..
include ../Makefile.internal

INC := -I../core_blas-qwrapper/ -I../control/ -I../runtime/ -I.. $(INC)

ALLSRC := cublasDtsmqr.cu cublasDpamm.cu cublasDparfb.cu cublasDlacpy.cu

ALLOBJ := $(ALLSRC:.cu=.o) 

all: $(LIBCOREBLASRW)

$(LIBCOREBLASRW): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean: 
	rm -f *.o

cleanall:
	rm -f *.o

#.cu.o :
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(INC) -c $< -o $@

.PHONY: all clean cleanall

