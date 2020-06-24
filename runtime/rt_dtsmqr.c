#include "runtime.h"
#include "core_blas-gpu.h"

//#pragma omp task inout([lda1*n1]A1, [lda2*n2]A2) in([ldt*nb]V, [ldt*nb]T) label(dtsmqr)
//CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);

#ifdef PLASMA_WITH_SMP
#pragma omp target device (smp) copy_deps
#pragma omp task inout([lda1*nb]A1, [lda2*nb]A2) in([ldv*nb]V, [ldt*nb]T) out([ib*nb]WORK) label(ldtsmqr_smp)
void CORE_ldtsmqr_ompss(PLASMA_enum side, PLASMA_enum trans, int m1, int n1, int m2, int n2, int k, int ib, int nb, double *A1, int lda1, double *A2, int lda2, const double *V, int ldv, const double *T, int ldt, double *WORK, int ldwork)
{
	CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}

#pragma omp target device (smp) copy_deps
#pragma omp task inout([lda1*nb]A1, [lda2*nb]A2) in([ldv*nb]V, [ldt*nb]T) out([ib*nb]WORK) label(rdtsmqr_smp)
void CORE_rdtsmqr_ompss(PLASMA_enum side, PLASMA_enum trans, int m1, int n1, int m2, int n2, int k, int ib, int nb, double *A1, int lda1, double *A2, int lda2, const double *V, int ldv, const double *T, int ldt, double *WORK, int ldwork)
{
	CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}
#endif


#ifdef PLASMA_WITH_CUDA_PURE
#pragma omp target device (cuda) copy_deps
#pragma omp task inout([lda1*nb]A1, [lda2*nb]A2) in([ldv*nb]V, [ldt*nb]T) out([ib*nb]WORK) label(ldtsmqr_cuda)
void CORE_ldtsmqr_ompss(PLASMA_enum side, PLASMA_enum trans, int m1, int n1, int m2, int n2, int k, int ib, int nb, double *A1, int lda1, double *A2, int lda2, const double *V, int ldv, const double *T, int ldt, double *WORK, int ldwork)
{
        /*
	CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
        */
	//printf("\n\n=============================> SALEM\n\n");
	cublasOperation_t cutrans;
	if ( trans == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDtsmqr(handle, CUBLAS_SIDE_LEFT, cutrans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}





#pragma omp target device (cuda) copy_deps
#pragma omp task inout([lda1*nb]A1, [lda2*nb]A2) in([ldv*nb]V, [ldt*nb]T) out([ib*nb]WORK) label(rdtsmqr_cuda)
void CORE_rdtsmqr_ompss(PLASMA_enum side, PLASMA_enum trans, int m1, int n1, int m2, int n2, int k, int ib, int nb, double *A1, int lda1, double *A2, int lda2, const double *V, int ldv, const double *T, int ldt, double *WORK, int ldwork)
{
/*
	CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
*/
	cublasOperation_t cutrans;
	if ( trans == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDtsmqr(handle, CUBLAS_SIDE_RIGHT, cutrans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}

#endif



#ifdef PLASMA_WITH_CUDA_HYBRID

#pragma omp target device (smp) copy_deps
#pragma omp task inout([lda1*nb]A1, [lda2*nb]A2) in([ldv*nb]V, [ldt*nb]T) out([ib*nb]WORK) label(ldtsmqr_hyb_smp)
void CORE_ldtsmqr_ompss(PLASMA_enum side, PLASMA_enum trans, int m1, int n1, int m2, int n2, int k, int ib, int nb, double *A1, int lda1, double *A2, int lda2, const double *V, int ldv, const double *T, int ldt, double *WORK, int ldwork)
{
	CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}

#pragma omp target device (smp) copy_deps
#pragma omp task inout([lda1*nb]A1, [lda2*nb]A2) in([ldv*nb]V, [ldt*nb]T) out([ib*nb]WORK) label(rdtsmqr_hyb_smp)
void CORE_rdtsmqr_ompss(PLASMA_enum side, PLASMA_enum trans, int m1, int n1, int m2, int n2, int k, int ib, int nb, double *A1, int lda1, double *A2, int lda2, const double *V, int ldv, const double *T, int ldt, double *WORK, int ldwork)
{
	CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}

//Alternative implementations

#pragma omp target device (cuda) copy_deps implements(CORE_ldtsmqr_ompss)
#pragma omp task inout([lda1*nb]A1, [lda2*nb]A2) in([ldv*nb]V, [ldt*nb]T) out([ib*nb]WORK) label(ldtsmqr_hyb_cuda)
void CORE_ldtsmqr_cuda(PLASMA_enum side, PLASMA_enum trans, int m1, int n1, int m2, int n2, int k, int ib, int nb, double *A1, int lda1, double *A2, int lda2, const double *V, int ldv, const double *T, int ldt, double *WORK, int ldwork)
{
        /*
	CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
        */
	//printf("\n\n=============================> SALEM\n\n");
	cublasOperation_t cutrans;
	if ( trans == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDtsmqr(handle, CUBLAS_SIDE_LEFT, cutrans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}


#pragma omp target device (cuda) copy_deps implements(CORE_rdtsmqr_ompss)
#pragma omp task inout([lda1*nb]A1, [lda2*nb]A2) in([ldv*nb]V, [ldt*nb]T) out([ib*nb]WORK) label(rdtsmqr_hyb_cuda)
void CORE_rdtsmqr_cuda(PLASMA_enum side, PLASMA_enum trans, int m1, int n1, int m2, int n2, int k, int ib, int nb, double *A1, int lda1, double *A2, int lda2, const double *V, int ldv, const double *T, int ldt, double *WORK, int ldwork)
{
/*
	CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
*/
	cublasOperation_t cutrans;
	if ( trans == PlasmaNoTrans )
		cutrans = CUBLAS_OP_N;
	else
		cutrans = CUBLAS_OP_T;

	cublasHandle_t handle = nanos_get_cublas_handle();
	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(handle, stream);
	cublasDtsmqr(handle, CUBLAS_SIDE_RIGHT, cutrans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}

#endif

void RT_CORE_dtsmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       double *A1, int lda1,
                       double *A2, int lda2,
                       const double *V, int ldv,
                       const double *T, int ldt)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dtsmqr(
			quark, task_flags,
			side, trans,
			m1, n1, m2, n2, k, ib, nb,
			A1, lda1,
			A2, lda2,
			V, ldv,
			T, ldt);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                /*
                */
		double *WORK = malloc(ib*nb*sizeof(double));
                #pragma omp register ([ib*nb]WORK)
                int ldwork = side == PlasmaLeft?ib:nb;
                
		//#pragma omp task inout([nb*nb]A1, [nb*nb]A2) in([nb*nb]V, [ib*nb]T) label(dtsmqr)
		//CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
		if (side == PlasmaLeft){
		   int ldwork = ib;
                   //printf("\n\n============> DTSMQR_L BEFORE m1 %d n1 %d m2 %d n2 %d k %d ib %d nb %d lda1 %d lda2 %d ldv %d ldt %d ldwork %d \n", m1, n1, m2, n2, k, ib, nb, lda1, lda2, ldv, ldt, ldwork);
		   CORE_ldtsmqr_ompss(side, trans, m1, n1, m2, n2, k, ib, nb, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);

		   //#pragma omp task inout([nb*nb]A1, [nb*nb]A2) in([nb*nb]V, [ib*nb]T) label(dtsmqr)
                   //#pragma omp task inout([lda1*nb]A1, [lda2*nb]A2) in([ldv*nb]V, [ldt*nb]T) label(ldtsmqr)
                   //CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
                   //printf("\n\n============> DTSMQR_L AFTER \n");
		} else {
		   int ldwork = nb;
                   //printf("\n\n============> DTSMQR_R BEFORE m1 %d n1 %d m2 %d n2 %d k %d ib %d nb %d lda1 %d lda2 %d ldv %d ldt %d ldwork %d \n", m1, n1, m2, n2, k, ib, nb, lda1, lda2, ldv, ldt, ldwork);
		   CORE_rdtsmqr_ompss(side, trans, m1, n1, m2, n2, k, ib, nb, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);

		   //#pragma omp task inout([nb*nb]A1, [nb*nb]A2) in([nb*nb]V, [ib*nb]T) label(dtsmqr)
                   //#pragma omp task inout([lda1*nb]A1, [lda2*nb]A2) in([ldv*nb]V, [ldt*nb]T) label(rdtsmqr)
                   //CORE_dtsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
                   //printf("\n\n============> DTSMQR_R AFTER \n");
		}
	}
}
