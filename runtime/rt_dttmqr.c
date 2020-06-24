#include "runtime.h"

void RT_CORE_dttmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       const double *V, int ldv,
                       const double *T, int ldt)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dttmqr(
			quark, task_flags,
			side, trans,
			m1, n1, m2, n2,
			k, ib, nb,
			A1, lda1,
			A2, lda2,
			V, ldv,
			T, ldt); 
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		double *WORK = malloc(ib*nb*sizeof(double));
		int ldwork;
                #pragma omp target device (smp) copy_deps
		#pragma omp task inout([nb*nb]A1, [nb*nb]A2) in([nb*nb]V, [ib*nb]T) label (dttmqr)
		CORE_dttmqr(side, trans, m1, n1, m2, n2, k, ib, 
				A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
	}
}
