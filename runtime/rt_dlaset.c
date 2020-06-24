#include "runtime.h"

void RT_CORE_dlaset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int M, int N,
                       double alpha, double beta,
                       double *A, int LDA)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
           QUARK_CORE_dlaset(
               quark, task_flags,
               uplo, M, N, alpha, beta,
               A, LDA);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task out([LDA*N]A) label(dlaset)
		LAPACKE_dlaset_work(
			LAPACK_COL_MAJOR,
			lapack_const(uplo),
			M, N, alpha, beta, A, LDA);
	}
}
