#include "runtime.h"

#ifdef PLASMA_WITH_SMP
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*k]A) inout([ldc*n]C) label(dsyrk_n_smp)
void CORE_dsyrk_ompss_n(PLASMA_enum uplo, PLASMA_enum trans, int n, int k, double alpha, double *A, int lda, double beta, double *C, int ldc, int nb)
{
	CORE_dsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([n*lda]A) inout([ldc*n]C) label(dsyrk_t_smp)
void CORE_dsyrk_ompss_t(PLASMA_enum uplo, PLASMA_enum trans, int n, int k, double alpha, double *A, int lda, double beta, double *C, int ldc, int nb)
{
	CORE_dsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
#endif


#ifdef PLASMA_WITH_CUDA_HYBRID
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*k]A) inout([ldc*n]C) label(dsyrk_n_hyb_smp)
void CORE_dsyrk_ompss_n(PLASMA_enum uplo, PLASMA_enum trans, int n, int k, double alpha, double *A, int lda, double beta, double *C, int ldc, int nb)
{
	CORE_dsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([n*lda]A) inout([ldc*n]C) label(dsyrk_t_hyb_smp)
void CORE_dsyrk_ompss_t(PLASMA_enum uplo, PLASMA_enum trans, int n, int k, double alpha, double *A, int lda, double beta, double *C, int ldc, int nb)
{
	CORE_dsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}

//Alternative implementations
#pragma omp target device (cuda) copy_deps implements(CORE_dsyrk_ompss_n)
#pragma omp task in([lda*k]A) inout([ldc*n]C) label(dsyrk_n_hyb_cuda)
void CORE_dsyrk_cuda_n(PLASMA_enum uplo, PLASMA_enum trans, int n, int k, double alpha, double *A, int lda, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( trans == PlasmaNoTrans)
		cutrans = CUBLAS_OP_N;
	else 
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsyrk(handle, cuplo, cutrans, n, k, &alpha, A, lda, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_dsyrk_ompss_t)
#pragma omp task in([n*lda]A) inout([ldc*n]C) label(dsyrk_t_hyb_cuda)
void CORE_dsyrk_cuda_t(PLASMA_enum uplo, PLASMA_enum trans, int n, int k, double alpha, double *A, int lda, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( trans == PlasmaNoTrans)
		cutrans = CUBLAS_OP_N;
	else 
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsyrk(handle, cuplo, cutrans, n, k, &alpha, A, lda, &beta, C, ldc);
}
#endif

#ifdef PLASMA_WITH_CUDA_PURE
#pragma omp target device (cuda) copy_deps
#pragma omp task in([lda*k]A) inout([ldc*n]C) label(dsyrk_n_cuda)
void CORE_dsyrk_ompss_n(PLASMA_enum uplo, PLASMA_enum trans, int n, int k, double alpha, double *A, int lda, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( trans == PlasmaNoTrans)
		cutrans = CUBLAS_OP_N;
	else 
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsyrk(handle, cuplo, cutrans, n, k, &alpha, A, lda, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps
#pragma omp task in([n*lda]A) inout([ldc*n]C) label(dsyrk_t_cuda)
void CORE_dsyrk_ompss_t(PLASMA_enum uplo, PLASMA_enum trans, int n, int k, double alpha, double *A, int lda, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	cublasOperation_t cutrans;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;
	if ( trans == PlasmaNoTrans)
		cutrans = CUBLAS_OP_N;
	else 
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsyrk(handle, cuplo, cutrans, n, k, &alpha, A, lda, &beta, C, ldc);
}
#endif

void RT_CORE_dsyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      double alpha, double *A, int lda,
                      double beta, double *C, int ldc)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dsyrk( quark, task_flags,
			uplo, trans, n, k, nb,
			alpha, A, lda, beta, C, ldc);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		if ( trans == PlasmaNoTrans ) {
			CORE_dsyrk_ompss_n(uplo, trans, n, k, alpha, A, lda, beta, C, ldc, nb);
		} else {
			CORE_dsyrk_ompss_t(uplo, trans, n, k, alpha, A, lda, beta, C, ldc, nb);
		}
	}
}
