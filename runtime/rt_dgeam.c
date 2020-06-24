#include "runtime.h"

#ifdef PLASMA_WITH_SMP
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*n]A) inout([ldb*n]B) label(dgeam)
void CORE_dgeam_ompss(  PLASMA_enum transA, PLASMA_enum transB,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double beta, double *B, int ldb)
{
	CORE_dgeam(transA, transB, m, n, nb, alpha, A, lda, beta, B, ldb);
}
#endif
// CUDA support (hybrid)
#ifdef PLASMA_WITH_CUDA_HYBRID
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*n]A) inout([ldb*n]B) label(dgeam)
void CORE_dgeam_ompss(  PLASMA_enum transA, PLASMA_enum transB,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double beta, double *B, int ldb)
{
	CORE_dgeam(transA, transB, m, n, nb, alpha, A, lda, beta, B, ldb);
}

//Alternative implementations
#pragma omp target device (cuda) copy_deps 
#pragma omp task in([lda*n]A) inout([ldb*n]B) label(dgeam)
void CORE_dgeam_cuda(  PLASMA_enum transA, PLASMA_enum transB,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double beta, double *B, int ldb)
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
        cublasDgeam(handle, trans0, trans1, m, n, &alpha, A, lda, &beta, B, ldb, B, ldb);
}
#endif

// CUDA support (pure) 
#ifdef PLASMA_WITH_CUDA_PURE
#pragma omp target device (cuda) copy_deps 
#pragma omp task in([lda*n]A) inout([ldb*n]B) label(dgeam)
void CORE_dgeam_ompss(  PLASMA_enum transA, PLASMA_enum transB,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double beta, double *B, int ldb)
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
        cublasDgeam(handle, trans0, trans1, m, n, &alpha, A, lda, &beta, B, ldb, B, ldb);
}
#endif

void RT_CORE_dgeam(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum transA, PLASMA_enum transB,
                      int m, int n, int nb,
                      double alpha, double *A, int lda,
                      double beta, double *B, int ldb)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dgeam(quark, task_flags, 
				transA, transB,
				m, n, nb, alpha, A, lda,
				beta, B, ldb);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		CORE_dgeam_ompss(transA, transB, m, n, nb, alpha, A, lda, beta, B, ldb);
	}
}


