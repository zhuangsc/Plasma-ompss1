#include "runtime.h"

void RT_CORE_dtsqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       double *T, int ldt)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dtsqrt(
			quark, task_flags,
			m, n, ib, nb,
			A1, lda1,
			A2, lda2,
			T, ldt);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		double *TAU = malloc(nb*sizeof(double));
		double *WORK = malloc(ib*nb*sizeof(double));
                //#pragma omp register ([ib*nb]WORK)
                //#pragma omp register ([nb]TAU)
		//#pragma omp task concurrent([nb*n]A1) inout([m*n]A2) out([ib*nb]T) label(dtsqrt)

                //printf("\n\n DTSQRT BEFORE m %d n %d ib %d lda1 %d lda2 %d ldt %d\n", m, n, ib, lda1, lda2, ldt);
                /*
                */
                #pragma omp target device (smp) copy_deps
		#pragma omp task inout([lda1*nb]A1) inout([lda2*nb]A2) out([ldt*nb]T) label(dtsqrt)
		CORE_dtsqrt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
                //printf("\n\n DTSQRT AFTER\n");
	}
}
