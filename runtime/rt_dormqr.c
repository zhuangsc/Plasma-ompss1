#include "runtime.h"

void RT_CORE_dormqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const double *A, int lda,
                       const double *T, int ldt,
                       double *C, int ldc)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dormqr(
			quark, task_flags,
			side, trans,
			m, n, k, ib, nb,
			A, lda,
			T, ldt,
			C, ldc);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		double *WORK = malloc(ib*nb*sizeof(double));
                //#pragma omp register ([ib*nb]WORK)
		//#pragma omp task concurrent([m*nb]A) in([ib*nb]T) inout([m*n]C) label(dormqr)

                //printf("\n\n DORMQR BEFORE m %d n %d k %d ib %d lda %d ldt %d ldc %d \n", m, n, k, ib, lda, ldt, ldc);
                /*
                */
                #pragma omp target device (smp) copy_deps
		#pragma omp task in([lda*nb]A, [ldt*nb]T) inout([ldc*nb]C) label(dormqr)
		CORE_dormqr(side, trans, m, n, k, ib, A, lda, T, ldt, C, ldc, WORK, nb);
                //printf("\n\n DORMQR AFTER\n");
	}
}
