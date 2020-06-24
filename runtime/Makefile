PLASMA_DIR= ..
include ../Makefile.internal

INC := -I../core_blas-qwrapper/ -I../control/ -I../core_blas-gpu $(INC)

ALLSRC= rt_dasum.c rt_dgemm.c rt_dgeam.c rt_dgeqrf.c rt_dgessq.c rt_dgetrf.c rt_dlacpy.c rt_dlange.c rt_dlaset.c rt_dlaswp.c rt_dormqr.c \
		rt_dplssq.c rt_dpotrf.c rt_dsyrk.c rt_dtrsm.c rt_dtsmqr.c rt_dtsqrt.c rt_dttmqr.c rt_dtrmm.c rt_dsyr2k.c rt_dsymm.c  rt_dlansy.c \
		rt_dtrtri.c rt_dlauum.c rt_dplgsy.c rt_dpltmg.c rt_dsyssq.c rt_dsrtdg.c rt_dswptr.c rt_dtwonm.c\
		runtime.c 

ALLOBJ = $(ALLSRC:.c=.o) 

all: $(LIBRUNTIME)

$(LIBRUNTIME): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean: 
	rm -f *.o

cleanall:
	rm -f *.o

.c.o :
	$(MCC) $(CFLAGS) $(INC) -c $< -o $@

.PHONY: all clean cleanall

