#include "runtime.h"

void RT_CORE_dtrtri(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum diag,
                       int n, int nb,
                       double *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dtrtri( quark, task_flags, uplo, diag, n, nb, A, lda, sequence, request, iinfo);
	} else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task inout([lda*n]A) label(dtrtri)
		CORE_dtrtri_rt(uplo, diag, n, A, lda);
	}
}

void CORE_dtrtri_rt(PLASMA_enum uplo, PLASMA_enum diag, int n, double *A, int lda)
{
	int info;
	info = LAPACKE_dtrtri_work(LAPACK_COL_MAJOR, lapack_const(uplo), lapack_const(diag), n, A, lda);
	//TODO plasma_sequence_flush
}
