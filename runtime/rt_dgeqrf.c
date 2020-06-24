#include "runtime.h"

void RT_CORE_dgeqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A, int lda,
                       double *T, int ldt)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
        QUARK_CORE_dgeqrt(
            quark, task_flags,
            m, n, ib, nb,
            A, lda,
            T, ldt);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                /*
                */
		double *TAU = malloc(nb*sizeof(double));
		double *WORK = malloc(ib*nb*sizeof(double));
                //#pragma omp register ([ib*nb]WORK)
                //#pragma omp register ([nb]TAU)
                //printf("\n\n DGEQRF BEFORE m %d n %d ib %d lda %d ldt %d\n", m, n, ib, lda, ldt);

                /*
                */
                #pragma omp target device (smp) copy_deps
		#pragma omp task inout([lda*nb]A) out([ldt*nb]T) label(dgeqrt)
		CORE_dgeqrt(m, n, ib, A, lda, T, ldt, TAU, WORK);
                //printf("\n\n DGEQRF AFTER\n");
	}
}
