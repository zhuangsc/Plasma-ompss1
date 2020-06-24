#include "runtime.h"

void RT_CORE_dsyssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const double *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF )
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();

	if ( plasma->runtime == PLASMA_QUARK ) {
		QUARK_CORE_dsyssq_f1(quark, task_flags, 
				uplo, n, A, lda, 
				scale, sumsq, 
				fake, szeF, paramF);
	} else if ( plasma->runtime == PLASMA_OMPSS) {
		if ( (fake == scale) && (paramF & GATHERV) ) {
                        #pragma omp target device (smp) copy_deps
			#pragma omp task in([lda*n]A) inout([1]scale, [1]sumsq) label(dsyssq)
			CORE_dsyssq(uplo, n, A, lda, scale, sumsq);
		} else {
                        #pragma omp target device (smp) copy_deps
			#pragma omp task in([lda*n]A) inout([1]scale, [1]sumsq) concurrent([szeF]fake) label(dsyssq)
			CORE_dsyssq(uplo, n, A, lda, scale, sumsq);
		}
	}
}
