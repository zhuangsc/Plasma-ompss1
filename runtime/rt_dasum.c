#include "runtime.h"

void RT_CORE_dasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                       const double *A, int lda, int szeA,
                       double *work, int szeW)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dasum(quark, task_flags, storev, uplo, M, N,
				A, lda, szeA, work, szeW);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task in([szeA]A) inout([szeW]work) label(dasum)
                
		CORE_dasum(storev, uplo, M, N, A, lda, work);
	}
}

void RT_CORE_dasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                          const double *A, int lda, int szeA,
                          double *work, int szeW, double *fake, int szeF)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dasum_f1(
			quark, task_flags,
			storev, uplo, M, N,
			A, lda, szeA,
			work, szeW,
			fake, szeF);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		if (work == fake) {
                        #pragma omp target device (smp) no_copy_deps
			#pragma omp task in([szeA]A) inout([szeW]work) label(dasum_f1) 
			CORE_dasum(storev, uplo, M, N, A, lda, work);
		} else {
                        #pragma omp target device (smp) no_copy_deps
			#pragma omp task in([szeA]A) inout([1]work) out([szeF]fake) label(dasum_f1)
			CORE_dasum(storev, uplo, M, N, A, lda, work);
		}
	}
}
