#include "runtime.h"

void RT_CORE_dlauum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       double *A, int lda)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dlauum(quark, task_flags,
			uplo, n, nb, A, lda);
	} else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task inout([lda*n]A) label(dlauum)
		CORE_dlauum_rt(uplo, n, A, lda);
	}
}

void CORE_dlauum_rt(PLASMA_enum uplo, int n, double *A, int lda)
{
	LAPACKE_dlauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), n, A, lda);
}
