#include "runtime.h"

void RT_CORE_dgessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, const double *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF )
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dgessq_f1(
			quark, task_flags,
			m, n,
			A, lda,
			scale, 
			sumsq,
			fake, szeF, paramF);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		if ((fake == scale) && (paramF & GATHERV)) {
                        #pragma omp target device (smp) copy_deps
			#pragma omp task in([lda*n]A) inout([1]scale, [1]sumsq) label(dgessq)
			CORE_dgessq( m, n, A, lda, scale, sumsq );
		} else {
                        #pragma omp target device (smp) copy_deps
                        #pragma omp task in([lda*n]A) inout([1]scale, [1]sumsq) concurrent(fake[0:szeF-1]) label(dgessq)
			CORE_dgessq( m, n, A, lda, scale, sumsq );
		}
	}
}
