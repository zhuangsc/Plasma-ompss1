#include "runtime.h"
#include "core_blas-gpu.h"

#ifdef PLASMA_WITH_SMP
#pragma omp target device (smp) copy_deps
#pragma omp task in([lda*(n-1)+m]A) out([ldb*(n-1)+m]B) label(dlacpy_smp)
void CORE_dlacpy_ompss(PLASMA_enum uplo, int m, int n, double *A, int lda, double *B, int ldb)
{
	LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, lapack_const(uplo),
				m, n, A, lda, B, ldb);
}
#endif

/*
*/
#ifdef PLASMA_WITH_CUDA_PURE
#pragma omp target device (cuda) copy_deps
#pragma omp task in([lda*(n-1)+m]A) out([ldb*(n-1)+m]B) label(dlacpy_cuda)
void CORE_dlacpy_ompss(PLASMA_enum uplo, int m, int n, double *A, int lda, double *B, int ldb)
{
        cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
        //printf("=========> cublasDlacpy m %d n %d lda %d ldb %d sizeA %d sizeB %d\n\n", m, n, lda, ldb, lda*(n-1)+m, ldb*(n-1)+m);
        cublasDlacpy(stream, uplo, m, n, A, lda, B, ldb );
	//LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, lapack_const(uplo),
				//m, n, A, lda, B, ldb);
}
#endif

void RT_CORE_dlacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int nb,
                       const double *A, int lda,
                       double *B, int ldb)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dlacpy(
			quark, task_flags, uplo,
			m, n, nb,
			A, lda,
			B, ldb);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		CORE_dlacpy_ompss(uplo, m, n, A, lda, B, ldb);
	}
}
