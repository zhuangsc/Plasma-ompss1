#include "runtime.h"

void RT_CORE_dpotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       double *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dpotrf(
			quark, task_flags,
			uplo, n, nb,
			A, lda,
			sequence, request, iinfo);

	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task inout([lda*n]A) label(dportf_smp)
		CORE_dpotrf(uplo, n, A, lda, &iinfo);
	}
}
