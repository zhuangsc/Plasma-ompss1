#include "runtime.h"

void RT_CORE_dlansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, PLASMA_enum uplo, int N,
                       const double *A, int LDA, int szeA,
                       int szeW, double *result)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dlansy( quark, task_flags, norm, uplo, N,
			A, LDA, szeA,
			szeW, result);
	} else if ( plasma->runtime == PLASMA_OMPSS ) {
		szeW = max(1, szeW);
		double *work = malloc(szeW * sizeof(double));
                #pragma omp target device (smp) copy_deps
		#pragma omp task in([szeA]A) out([1]result) label(dlansy)
		*result = LAPACKE_dlansy_work(LAPACK_COL_MAJOR, lapack_const(norm),
				lapack_const(uplo), N, A, LDA, work);
	}
}

void RT_CORE_dlansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, PLASMA_enum uplo, int N,
                          const double *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dlansy_f1( quark, task_flags, norm, uplo, N,
			A, LDA, szeA,
			szeW, result,
			fake, szeF);
	} else if ( plasma->runtime == PLASMA_OMPSS ) {
		szeW = max(1, szeW);
		double *work = malloc(szeW * sizeof(double));
		if ( result == fake ) {
                        #pragma omp target device (smp) copy_deps
			#pragma omp task in([szeA]A) out([1]result) label(dlansy_f1)
			*result = LAPACKE_dlansy_work(LAPACK_COL_MAJOR, lapack_const(norm),
					lapack_const(uplo), N, A, LDA, work);
		} else {
                        #pragma omp target device (smp) copy_deps
			#pragma omp task in([szeA]A) out([1]result) fake([szeF]fake) label(dlansy_f1)
			*result = LAPACKE_dlansy_work(LAPACK_COL_MAJOR, lapack_const(norm),
					lapack_const(uplo), N, A, LDA, work);
		}
	}
}
