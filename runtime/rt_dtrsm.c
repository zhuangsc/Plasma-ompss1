#include "runtime.h"

#ifdef PLASMA_WITH_SMP
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*m]A) inout([ldb*n]B) label(ldtrsm_smp)
void CORE_ldtrsm_ompss(PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int m, int n, double alpha, double *A, int lda, double *B, int ldb, int nb)
{
	CORE_dtrsm(PlasmaLeft, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*n]A) inout([ldb*n]B) label(rdtrsm_smp)
void CORE_rdtrsm_ompss(PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int m, int n, double alpha, double *A, int lda, double *B, int ldb, int nb)
{
	CORE_dtrsm(PlasmaRight, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
#endif


#ifdef PLASMA_WITH_CUDA_HYBRID
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*m]A) inout([ldb*n]B) label(ldtrsm_hyb_smp)
void CORE_ldtrsm_ompss(PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int m, int n, double alpha, double *A, int lda, double *B, int ldb, int nb)
{
	CORE_dtrsm(PlasmaLeft, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*n]A) inout([ldb*n]B) label(rdtrsm_hyb_smp)
void CORE_rdtrsm_ompss(PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int m, int n, double alpha, double *A, int lda, double *B, int ldb, int nb)
{
	CORE_dtrsm(PlasmaRight, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}

//Alternative implementations
#pragma omp target device (cuda) copy_deps implements(CORE_ldtrsm_ompss)
#pragma omp task in([lda*m]A) inout([ldb*n]B) label(ldtrsm_hyb_cuda)
void CORE_ldtrsm_cuda(PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int m, int n, double alpha, double *A, int lda, double *B, int ldb, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	cublasDiagType_t cudiag;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( transA == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;
	if ( diag == PlasmaNonUnit )
		cudiag = CUBLAS_DIAG_NON_UNIT;
	else
		cudiag = CUBLAS_DIAG_UNIT;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDtrsm(handle, CUBLAS_SIDE_LEFT, cuplo, cutrans, cudiag, m, n, &alpha, A, lda, B, ldb);
}

#pragma omp target device (cuda) copy_deps implements(CORE_rdtrsm_ompss)
#pragma omp task in([lda*n]A) inout([ldb*n]B) label(rdtrsm_hyb_cuda)
void CORE_rdtrsm_cuda(PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int m, int n, double alpha, double *A, int lda, double *B, int ldb, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	cublasDiagType_t cudiag;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( transA == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;
	if ( diag == PlasmaNonUnit )
		cudiag = CUBLAS_DIAG_NON_UNIT;
	else
		cudiag = CUBLAS_DIAG_UNIT;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDtrsm(handle, CUBLAS_SIDE_RIGHT, cuplo, cutrans, cudiag, m, n, &alpha, A, lda, B, ldb);
}
#endif


#ifdef PLASMA_WITH_CUDA_PURE
#pragma omp target device (cuda) copy_deps
#pragma omp task in([lda*m]A) inout([ldb*n]B) label(ldtrsm_cuda)
void CORE_ldtrsm_ompss(PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int m, int n, double alpha, double *A, int lda, double *B, int ldb, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	cublasDiagType_t cudiag;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( transA == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;
	if ( diag == PlasmaNonUnit )
		cudiag = CUBLAS_DIAG_NON_UNIT;
	else
		cudiag = CUBLAS_DIAG_UNIT;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDtrsm(handle, CUBLAS_SIDE_LEFT, cuplo, cutrans, cudiag, m, n, &alpha, A, lda, B, ldb);
}

#pragma omp target device (cuda) copy_deps 
#pragma omp task in([lda*n]A) inout([ldb*n]B) label(rdtrsm_cuda)
void CORE_rdtrsm_ompss(PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, int m, int n, double alpha, double *A, int lda, double *B, int ldb, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	cublasDiagType_t cudiag;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( transA == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;
	if ( diag == PlasmaNonUnit )
		cudiag = CUBLAS_DIAG_NON_UNIT;
	else
		cudiag = CUBLAS_DIAG_UNIT;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDtrsm(handle, CUBLAS_SIDE_RIGHT, cuplo, cutrans, cudiag, m, n, &alpha, A, lda, B, ldb);
}
#endif

void RT_CORE_dtrsm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dtrsm(
			quark, task_flags,
			side, uplo, transA, diag,
			m, n, nb,
			alpha, A, lda, B, ldb);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
//		cublasSideMode_t cuside;
		if (side == PlasmaLeft){
			CORE_ldtrsm_ompss(uplo, transA, diag, m, n, alpha, A, lda, B, ldb, nb);
		} else {
			CORE_rdtrsm_ompss(uplo, transA, diag, m, n, alpha, A, lda, B, ldb, nb);
		}
	}
}
