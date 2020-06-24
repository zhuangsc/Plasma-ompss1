#include "runtime.h"

void RT_CORE_dplgsy( Quark *quark, Quark_Task_Flags *task_flags,
                        double bump, int m, int n, double *A, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed )
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dplgsy(quark, task_flags,
			bump, m, n, A, lda,
			bigM, m0, n0, seed );
	} else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task out([lda*n]A) label(dplgsy)
		CORE_dplgsy_rt(bump, m, n, A, lda, bigM, m0, n0, seed);
	}
}

void CORE_dplgsy_rt(int bump, int m, int n, double *A, int lda, int bigM,
		int m0, int n0, unsigned long long int seed)
{
	CORE_dplgsy(bump, m, n, A, lda, bigM, m0, n0, seed);
}
