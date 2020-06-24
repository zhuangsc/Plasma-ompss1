#include "runtime.h"

#ifdef PLASMA_WITH_SMP
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*k]A, [ldb*n]B) inout([ldc*n]C) label(dgemm_nn_smp)
void CORE_dgemm_ompss_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*k]A, [k*ldb]B) inout([ldc*n]C) label(dgemm_nt_smp)
void CORE_dgemm_ompss_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([m*lda]A, [ldb*n]B) inout([ldc*n]C) label(dgemm_tn_smp)
void CORE_dgemm_ompss_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([m*lda]A, [k*ldb]B) inout([ldc*n]C) label(dgemm_tt_smp)
void CORE_dgemm_ompss_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}
#endif

/*-------------------------------------------------------------------------------------*/
// CUDA support (hybrid)
#ifdef PLASMA_WITH_CUDA_HYBRID
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*k]A, [ldb*n]B) inout([ldc*n]C) label(dgemm_nn_smp)
void CORE_dgemm_ompss_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp)  copy_deps
#pragma omp task in([lda*k]A, [k*ldb]B) inout([ldc*n]C) label(dgemm_nt_smp)
void CORE_dgemm_ompss_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([m*lda]A, [ldb*n]B) inout([ldc*n]C) label(dgemm_tn_smp)
void CORE_dgemm_ompss_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) copy_deps 
#pragma omp task in([m*lda]A, [k*ldb]B) inout([ldc*n]C) label(dgemm_tt_smp)
void CORE_dgemm_ompss_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

//Alternative implementations
#pragma omp target device (cuda) copy_deps implements(CORE_dgemm_ompss_nn)
#pragma omp task in([lda*k]A, [ldb*n]B) inout([ldc*n]C) label(dgemm_nn_cuda)
void CORE_dgemm_cuda_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_dgemm_ompss_nt)
#pragma omp task in([lda*k]A, [k*ldb]B) inout([ldc*n]C) label(dgemm_nt_cuda)
void CORE_dgemm_cuda_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_dgemm_ompss_tn)
#pragma omp task in([m*lda]A, [ldb*n]B) inout([ldc*n]C) label(dgemm_tn_cuda)
void CORE_dgemm_cuda_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_dgemm_ompss_tt)
#pragma omp task in([m*lda]A, [k*ldb]B) inout([ldc*n]C) label(dgemm_tt_cuda)
void CORE_dgemm_cuda_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}
#endif

/*-------------------------------------------------------------------------------------*/
// CUDA support (pure) 
#ifdef PLASMA_WITH_CUDA_PURE
#pragma omp target device (cuda) copy_deps 
#pragma omp task in([lda*k]A, [ldb*n]B) inout([ldc*n]C) label(dgemm_nn_cuda)
void CORE_dgemm_ompss_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps
#pragma omp task in([lda*k]A, [k*ldb]B) inout([ldc*n]C) label(dgemm_nt_cuda)
void CORE_dgemm_ompss_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps
#pragma omp task in([m*lda]A, [ldb*n]B) inout([ldc*n]C) label(dgemm_tn_cuda)
void CORE_dgemm_ompss_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps
#pragma omp task in([m*lda]A, [k*ldb]B) inout([ldc*n]C) label(dgemm_tt_cuda)
void CORE_dgemm_ompss_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}
#endif

//-------------------------------------------------
//concurrent dgemm2 tasks
//-------------------------------------------------
//#ifdef PLASMA_WITH_SMP
#pragma omp target device (smp) no_copy_deps
#pragma omp task in([lda*k]A, [ldb*n]B) concurrent([ldc*n]C) label(dgemm2_nn_smp)
void CORE_dgemm2_ompss_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) no_copy_deps
#pragma omp task in([lda*k]A, [k*ldb]B) concurrent([ldc*n]C) label(dgemm2_nt_smp)
void CORE_dgemm2_ompss_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) no_copy_deps
#pragma omp task in([m*lda]A, [ldb*n]B) concurrent([ldc*n]C) label(dgemm2_tn_smp)
void CORE_dgemm2_ompss_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) no_copy_deps
#pragma omp task in([m*lda]A, [k*ldb]B) concurrent([ldc*n]C) label(dgemm2_tt_smp)
void CORE_dgemm2_ompss_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}
//#endif

/*-------------------------------------------------------------------------------------*/
// CUDA support (hybrid)
/*
#ifdef PLASMA_WITH_CUDA_HYBRID
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*k]A, [ldb*n]B) concurrent([ldc*n]C) label(dgemm2_nn_smp)
void CORE_dgemm2_ompss_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp)  copy_deps
#pragma omp task in([lda*k]A, [k*ldb]B) concurrent([ldc*n]C) label(dgemm2_nt_smp)
void CORE_dgemm2_ompss_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([m*lda]A, [ldb*n]B) concurrent([ldc*n]C) label(dgemm2_tn_smp)
void CORE_dgemm2_ompss_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) copy_deps 
#pragma omp task in([m*lda]A, [k*ldb]B) concurrent([ldc*n]C) label(dgemm2_tt_smp)
void CORE_dgemm2_ompss_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

//Alternative implementations
#pragma omp target device (cuda) copy_deps implements(CORE_dgemm2_ompss_nn)
#pragma omp task in([lda*k]A, [ldb*n]B) concurrent([ldc*n]C) label(dgemm2_nn_cuda)
void CORE_dgemm2_cuda_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_dgemm2_ompss_nt)
#pragma omp task in([lda*k]A, [k*ldb]B) concurrent([ldc*n]C) label(dgemm2_nt_cuda)
void CORE_dgemm2_cuda_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_dgemm2_ompss_tn)
#pragma omp task in([m*lda]A, [ldb*n]B) concurrent([ldc*n]C) label(dgemm2_tn_cuda)
void CORE_dgemm2_cuda_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_dgemm2_ompss_tt)
#pragma omp task in([m*lda]A, [k*ldb]B) concurrent([ldc*n]C) label(dgemm2_tt_cuda)
void CORE_dgemm2_cuda_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}
#endif
*/

/*-------------------------------------------------------------------------------------*/
// CUDA support (pure) 
/*
#ifdef PLASMA_WITH_CUDA_PURE
#pragma omp target device (cuda) copy_deps 
#pragma omp task in([lda*k]A, [ldb*n]B) concurrent([ldc*n]C) label(dgemm2_nn_cuda)
void CORE_dgemm2_ompss_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps
#pragma omp task in([lda*k]A, [k*ldb]B) concurrent([ldc*n]C) label(dgemm2_nt_cuda)
void CORE_dgemm2_ompss_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps
#pragma omp task in([m*lda]A, [ldb*n]B) concurrent([ldc*n]C) label(dgemm2_tn_cuda)
void CORE_dgemm2_ompss_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps
#pragma omp task in([m*lda]A, [k*ldb]B) concurrent([ldc*n]C) label(dgemm2_tt_cuda)
void CORE_dgemm2_ompss_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}
#endif
*/


//----------------------------------------------------------
//dgemm_f2 tasks
//----------------------------------------------------------
//#ifdef PLASMA_WITH_SMP
#pragma omp target device (smp) no_copy_deps
#pragma omp task in([lda*k]A, [ldb*n]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_nn_smp)
void CORE_dgemm_f2_ompss_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) no_copy_deps
#pragma omp task in([lda*k]A, [k*ldb]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_nt_smp)
void CORE_dgemm_f2_ompss_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) no_copy_deps
#pragma omp task in([m*lda]A, [ldb*n]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_tn_smp)
void CORE_dgemm_f2_ompss_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) no_copy_deps
#pragma omp task in([m*lda]A, [k*ldb]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_tt_smp)
void CORE_dgemm_f2_ompss_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}
//#endif

/*-------------------------------------------------------------------------------------*/
// CUDA support (hybrid)
/*
#ifdef PLASMA_WITH_CUDA_HYBRID
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*k]A, [ldb*n]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_nn_smp)
void CORE_dgemm_f2_ompss_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*k]A, [k*ldb]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_nt_smp)
void CORE_dgemm_f2_ompss_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([m*lda]A, [ldb*n]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_tn_smp)
void CORE_dgemm_f2_ompss_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

#pragma omp target device (smp) copy_deps
#pragma omp task in([m*lda]A, [k*ldb]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_tt_smp)
void CORE_dgemm_f2_ompss_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	CORE_dgemm( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc);
}

//Alternative implementations
#pragma omp target device (cuda) copy_deps implements(CORE_dgemm_f2_ompss_nn)
#pragma omp task in([lda*k]A, [ldb*n]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_nn_cuda)
void CORE_dgemm_f2_cuda_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_dgemm_f2_ompss_nt)
#pragma omp task in([lda*k]A, [k*ldb]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_nt_cuda)
void CORE_dgemm_f2_cuda_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_dgemm_f2_ompss_tn)
#pragma omp task in([m*lda]A, [ldb*n]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_tn_cuda)
void CORE_dgemm_f2_cuda_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps implements(CORE_dgemm_f2_ompss_tt)
#pragma omp task in([m*lda]A, [k*ldb]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_tt_cuda)
void CORE_dgemm_f2_cuda_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}
#endif
*/

/*-------------------------------------------------------------------------------------*/
// CUDA support (pure) 
/*
#ifdef PLASMA_WITH_CUDA_PURE
#pragma omp target device (cuda) copy_deps 
#pragma omp task in([lda*k]A, [ldb*n]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_nn_cuda)
void CORE_dgemm_f2_ompss_nn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps 
#pragma omp task in([lda*k]A, [k*ldb]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_nt_cuda)
void CORE_dgemm_f2_ompss_nt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps 
#pragma omp task in([m*lda]A, [ldb*n]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_tn_cuda)
void CORE_dgemm_f2_ompss_tn(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#pragma omp target device (cuda) copy_deps 
#pragma omp task in([m*lda]A, [k*ldb]B) inout([ldc*n]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2_tt_cuda)
void CORE_dgemm_f2_ompss_tt(PLASMA_enum transA, PLASMA_enum transB, int m , int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc, int nb, double *fake1, int szefake1, double *fake2, int szefake2)
{
	cublasOperation_t trans0, trans1;
	if ( transA == PlasmaNoTrans)
		trans0 = CUBLAS_OP_N;
	else
		trans0 = CUBLAS_OP_T;
	if ( transB == PlasmaNoTrans)
		trans1 = CUBLAS_OP_N;
	else
		trans1 = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDgemm(handle, trans0, trans1, m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}

#endif
*/


void RT_CORE_dgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum transA, int transB,
                      int m, int n, int k, int nb,
                      double alpha, double *A, int lda,
					  double *B, int ldb,
                      double beta, double *C, int ldc)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dgemm(quark, task_flags, 
				transA, transB,
				m, n, k, nb, alpha, A, lda, B, ldb,
				beta, C, ldc);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		if (transA == PlasmaNoTrans){
			if (transB == PlasmaNoTrans){
				CORE_dgemm_ompss_nn( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb);
			} else {
				CORE_dgemm_ompss_nt( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb);
			}
		} else {
			if (transB == PlasmaNoTrans){
				CORE_dgemm_ompss_tn( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb);
			} else {
				CORE_dgemm_ompss_tt( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb);
			}
		}
	}
}

void RT_CORE_dgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum transA, int transB,
                        int m, int n, int k, int nb,
                        double alpha, const double *A, int lda,
                        const double *B, int ldb,
                        double beta, double *C, int ldc)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dgemm2(
			quark, task_flags,
			transA, transB,
			m, n, k, nb,
			alpha, A, lda,
			B, ldb,
			beta,  C, ldc);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		if (transA == PlasmaNoTrans){
			if (transB == PlasmaNoTrans){
				CORE_dgemm2_ompss_nn( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb);
			} else {
				CORE_dgemm2_ompss_nt( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb);
			}
		} else {
			if (transB == PlasmaNoTrans){
				CORE_dgemm2_ompss_tn( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb);
			} else {
				CORE_dgemm2_ompss_tt( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb);
			}
		}
	}
}

void RT_CORE_dgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, int transB,
                         int m, int n, int k, int nb,
                         double alpha, const double *A, int lda,
					   	 const double *B, int ldb,
                         double beta, double *C, int ldc,
                         double *fake1, int szefake1, int flag1,
                         double *fake2, int szefake2, int flag2)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dgemm_f2(
			quark, task_flags,
			transA, transB,
			m, n, k, nb,
			alpha, A, lda,
				   B, ldb,
			beta,  C, ldc,
			fake1, szefake1, flag1,
			fake2, szefake2, flag2);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		if (transA == PlasmaNoTrans){
			if (transB == PlasmaNoTrans){
				CORE_dgemm_f2_ompss_nn( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb, fake1, szefake1, fake2, szefake2);
			} else {
				CORE_dgemm_f2_ompss_nt( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb, fake1, szefake1, fake2, szefake2);
			}
		} else {
			if (transB == PlasmaNoTrans){
				CORE_dgemm_f2_ompss_tn( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb, fake1, szefake1, fake2, szefake2);
			} else {
				CORE_dgemm_f2_ompss_tt( transA, transB, m, n, k, (alpha), A, lda, B, ldb, (beta), C, ldc, nb, fake1, szefake1, fake2, szefake2);
			}
		}
//		printf("dgemm_f2 A:%p, B: %p, C: %p, fake1: %p, fake2: %p\n", A, B, C, fake1, &fake2);
//		#pragma omp task in([1]A, [1]B) inout([1]C) inout([1]fake1) in([szefake2]fake2) label(dgemm_f2)
//		CORE_dgemm(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
	}
}
