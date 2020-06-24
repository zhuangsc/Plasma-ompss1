#include "runtime.h"

#ifdef PLASMA_WITH_SMP
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*m]A, [ldb*n]B) inout([ldc*n]C) label(ldsymm_smp)
void CORE_ldsymm_ompss(PLASMA_enum uplo, int m , int n, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dsymm(PlasmaLeft, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*n]A, [ldb*n]B) inout([ldc*n]C) label(rdsymm_smp)
void CORE_rdsymm_ompss(PLASMA_enum uplo, int m , int n, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dsymm(PlasmaRight, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
#endif


//CUDA support
#ifdef PLASMA_WITH_CUDA_HYBRID
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*m]A, [ldb*n]B) inout([ldc*n]C) label(ldsymm_hyb_smp)
void CORE_ldsymm_ompss(PLASMA_enum uplo, int m , int n, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dsymm(PlasmaLeft, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*n]A, [ldb*n]B) inout([ldc*n]C) label(rdsymm_hyb_smp)
void CORE_rdsymm_ompss(PLASMA_enum uplo, int m , int n, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dsymm(PlasmaRight, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}

//Alternative implementations
#pragma omp target device (cuda) copy_deps implements(CORE_ldsymm_ompss)
#pragma omp task in([lda*m]A, [ldb*n]B) inout([ldc*n]C) label(ldsymm_hyb_cuda)
void CORE_ldsymm_cuda(PLASMA_enum uplo, int m , int n, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsymm(handle, CUBLAS_SIDE_LEFT, cuplo, m, n, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_rdsymm_ompss)
#pragma omp task in([lda*n]A, [ldb*n]B) inout([ldc*n]C) label(rdsymm_hyb_cuda)
void CORE_rdsymm_cuda(PLASMA_enum uplo, int m , int n, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsymm(handle, CUBLAS_SIDE_RIGHT, cuplo, m, n, &alpha, A, lda, B, ldb, &beta, C, ldc);
}
#endif

#ifdef PLASMA_WITH_CUDA_PURE
#pragma omp target device (cuda) copy_deps
#pragma omp task in([lda*m]A, [ldb*n]B) inout([ldc*n]C) label(ldsymm_cuda)
void CORE_ldsymm_ompss(PLASMA_enum uplo, int m , int n, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsymm(handle, CUBLAS_SIDE_LEFT, cuplo, m, n, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps
#pragma omp task in([lda*n]A, [ldb*n]B) inout([ldc*n]C) label(rdsymm_cuda)
void CORE_rdsymm_ompss(PLASMA_enum uplo, int m , int n, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasFillMode_t cuplo;
	if ( uplo == PlasmaUpper )
		cuplo = CUBLAS_FILL_MODE_UPPER;
	else
		cuplo = CUBLAS_FILL_MODE_LOWER;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDsymm(handle, CUBLAS_SIDE_RIGHT, cuplo, m, n, &alpha, A, lda, B, ldb, &beta, C, ldc);
}
#endif

void RT_CORE_dsymm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double *B, int ldb,
                      double beta, double *C, int ldc)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dsymm( quark, task_flags,
			side, uplo, m, n, nb,
			alpha, A, lda, B, ldb,
			beta, C, ldc);
	} else if (plasma->runtime == PLASMA_OMPSS) {
		if (side == PlasmaLeft){
			CORE_ldsymm_ompss(uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc, nb);
		} else {
			CORE_rdsymm_ompss(uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc, nb);
		}
	}
}
