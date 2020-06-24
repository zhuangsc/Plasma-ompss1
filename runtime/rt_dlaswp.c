#include "runtime.h"

void RT_CORE_dlaswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, double *A, int lda,
                       int i1,  int i2, const int *ipiv, int inc)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dlaswp(
			quark, task_flags,
			n, A, lda, i1, i2, ipiv, inc);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task inout([lda*n]A) in([n]ipiv) label(dlaswp)
		LAPACKE_dlaswp_work(LAPACK_COL_MAJOR, n, A, lda, i1, i2, ipiv, inc);
	}
}

void RT_CORE_dlaswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, double *A, int lda,
                          int i1,  int i2, const int *ipiv, int inc,
                          double *fake1, int szefake1, int flag1,
                          double *fake2, int szefake2, int flag2)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dlaswp_f2(
		quark, task_flags,
		n, A, lda, i1, i2, ipiv, inc,
		fake1, szefake1, flag1,
		fake2, szefake2, flag2 );
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
                #pragma omp task inout([lda*n]A, [szefake2]fake2) in([n]ipiv, [szefake1]fake1) label(dlaswp_f2)
		LAPACKE_dlaswp_work(LAPACK_COL_MAJOR, n, A, lda, i1, i2, ipiv, inc);
	}
}

void RT_CORE_dlaswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, double *Aij,
                                 int i1,  int i2, const int *ipiv, int inc,
                                 double *fake1, int szefake1, int flag1,
                                 double *fake2, int szefake2, int flag2)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dlaswp_ontile_f2(
		quark, task_flags,
		descA, Aij, i1, i2, ipiv, inc,
		fake1, szefake1, flag1,
		fake2, szefake2, flag2 );
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
//		printf("dlaswp_ontile_f2 Aij: %p, fake1: %p, fake2: %p\n", Aij, fake1, fake2);
                #pragma omp target device (smp) no_copy_deps
		#pragma omp task in([(i2-i1+1)*abs(inc)]ipiv) inout([1]Aij) in([1]fake1) inout([szefake2]fake2) label(dlaswp_ontile_f2)
		CORE_dlaswp_ontile(descA, i1, i2, ipiv, inc);
	}
}
