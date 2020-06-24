#include "runtime.h"

void RT_CORE_dgetrf(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       double *A, int lda,
                       int *IPIV,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       PLASMA_bool check_info, int iinfo)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dgetrf(
			quark, task_flags,
			m, n, nb,
			A, lda, IPIV,
			sequence, request, check_info, iinfo );
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task inout([nb*nb]A) out([nb]IPIV) label(dgetrf)
		RT_LAPACKE_dgetrf_work(m, n, A, lda, IPIV);
	}
}

void RT_LAPACKE_dgetrf_work(int m, int n, double *A, int lda, int *IPIV)
{
	int info;
	info = LAPACKE_dgetrf_work(LAPACK_COL_MAJOR, m, n, A, lda, IPIV);
	if (info != PLASMA_SUCCESS) {
		int i;
		for(i=info-1; i<min(m,n); i++)
			IPIV[i] = i+1;
	}
}

void RT_CORE_dgetrf_reclap(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              double *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dgetrf_reclap(
			quark, task_flags,
			m, n, nb,
			A, lda, IPIV,
			sequence, request, check_info, iinfo,
			nbthread );
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task inout([nb*nb]A) out([nb]IPIV) label(dgetrf_reclap)
		CORE_dgetrf_reclap(m, n, A, lda, IPIV, iinfo);
	}
}

void RT_CORE_dgetrf_rectil(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A, double *Amn, int size,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dgetrf_rectil(
			quark, task_flags,
			A, Amn, size, IPIV,
			sequence, request, check_info, iinfo,
			nbthread );
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) no_copy_deps
		#pragma omp task inout([1]Amn) out([A.n]IPIV) label(dgetrf_rectil)
		RT_CORE_dgetrf_rectil_info(A, IPIV, nbthread);
	}
}

void RT_CORE_dgetrf_rectil_info(const PLASMA_desc A, int *IPIV, int nbthread)
{
	int info[3];
	info[1] = 0;
//	info[1] = RT_local_rank(t_id);
	info[2] = 1;
	CORE_dgetrf_rectil(A, IPIV, info);
}
