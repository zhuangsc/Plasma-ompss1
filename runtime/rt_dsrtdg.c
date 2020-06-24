#include "runtime.h"

void RT_dsrtdg(Quark *quark, Quark_Task_Flags *task_flags,
                       int M, const double *A, int lda,
                       double *work, int ldw)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();

        #pragma omp target device (smp) copy_deps
        #pragma omp task  in([M*lda]A) out([M]work) label(dsrtdg)
        RT_dsrtdg2( M, A, lda, work, ldw);
}


void RT_dsrtdg2( int M, const double *A, int lda,
               double *work, int ldw)
{
        int i, j;
        {
           for(i = 0; i < M; i++){
             work[i] = A[i*lda+i];
           }
        }
}


void RT_sort( int M, double *work)
{        
    #pragma omp target device (smp) copy_deps
    #pragma omp task inout([M]work) label(sort dsrtdg)
    LAPACKE_dlasrt_work( 'I', M, work);

}
