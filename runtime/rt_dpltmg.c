#include "runtime.h"

void RT_CORE_dpltmg( Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum mtxtype, int m, int n, double *A, int lda,
                         int gM, int gN, int m0, int n0, unsigned long long int seed )
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dpltmg(quark, task_flags,
			mtxtype, m, n, A, lda,
			gM, gN, m0, n0, seed );
	} else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task out([lda*n]A) label(dpltmg)
		CORE_dpltmg_rt(mtxtype, m, n, A, lda, gM, gN, m0, n0, seed);
	}
}

void CORE_dpltmg_rt(int mtxtype, int m, int n, double *A, int lda, int gM, int gN, int m0, int n0, unsigned long long int seed)
{
	CORE_dpltmg(mtxtype, m, n, A, lda, gM, gN, m0, n0, seed);
}
