#include "runtime.h"

#ifdef PLASMA_WITH_SMP
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*k]A, [ldb*k]B) inout([ldc*n]C) label(dsyr2k_n_smp)
void CORE_dsyr2k_ompss_n(PLASMA_enum uplo, PLASMA_enum trans, int n , int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([n*lda]A, [n*ldb]B) inout([ldc*n]C) label(dsyr2k_t_smp)
void CORE_dsyr2k_ompss_t(PLASMA_enum uplo, PLASMA_enum trans, int n , int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
#endif


#ifdef PLASMA_WITH_CUDA_HYBRID
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*k]A, [ldb*k]B) inout([ldc*n]C) label(dsyr2k_n_hyb_smp)
void CORE_dsyr2k_ompss_n(PLASMA_enum uplo, PLASMA_enum trans, int n , int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([n*lda]A, [n*ldb]B) inout([ldc*n]C) label(dsyr2k_t_hyb_smp)
void CORE_dsyr2k_ompss_t(PLASMA_enum uplo, PLASMA_enum trans, int n , int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

//Alternative implementations
#pragma omp target device (cuda) copy_deps implements(CORE_dsyr2k_ompss_n)
#pragma omp task in([lda*k]A, [ldb*k]B) inout([ldc*n]C) label(dsyr2k_n_hyb_cuda)
void CORE_dsyr2k_cuda_n(PLASMA_enum uplo, PLASMA_enum trans, int n , int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( trans == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsyr2k(handle, cuplo, cutrans, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_dsyr2k_ompss_t)
#pragma omp task in([n*lda]A, [n*ldb]B) inout([ldc*n]C) label(dsyr2k_t_hyb_cuda)
void CORE_dsyr2k_cuda_t(PLASMA_enum uplo, PLASMA_enum trans, int n , int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( trans == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsyr2k(handle, cuplo, cutrans, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}
#endif

#ifdef PLASMA_WITH_CUDA_PURE
#pragma omp target device (cuda) copy_deps 
#pragma omp task in([lda*k]A, [ldb*k]B) inout([ldc*n]C) label(dsyr2k_n_cuda)
void CORE_dsyr2k_ompss_n(PLASMA_enum uplo, PLASMA_enum trans, int n , int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( trans == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsyr2k(handle, cuplo, cutrans, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps
#pragma omp task in([n*lda]A, [n*ldb]B) inout([ldc*n]C) label(dsyr2k_t_cuda)
void CORE_dsyr2k_ompss_t(PLASMA_enum uplo, PLASMA_enum trans, int n , int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( trans == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsyr2k(handle, cuplo, cutrans, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}
#endif

void RT_CORE_dsyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       double alpha, double *A, int lda,
                       double *B, int ldb,
                       double beta, double *C, int ldc)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dsyr2k( quark, task_flags,
			uplo, trans, n, k, nb,
			alpha, A, lda, B, ldb,
			beta, C, ldc);
	} else if (plasma->runtime == PLASMA_OMPSS) {
		if (trans == PlasmaNoTrans) {
			CORE_dsyr2k_ompss_n(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc, nb);
		} else {
			CORE_dsyr2k_ompss_t(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc, nb);
		}
	}
}
