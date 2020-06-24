#include "runtime.h"

void CORE_dplssq(int m, double *A, double *result)
{
	int i;
    for( i=1; i<m; i++ ) {
        if( A[0] < A[i] ) {
            A[1] = A[2*i+1] + (A[1]  * (( A[0] / A[2*i] ) * ( A[0] / A[2*i] )));
            A[0] = A[2*i];
        } else {
            A[1] = A[1]  + (A[2*i+1] * (( A[2*i] / A[0] ) * ( A[2*i] / A[0] )));
        }
    }

    *result = A[0] * sqrt( A[1] );
}

void RT_CORE_dplssq( Quark *quark, Quark_Task_Flags *task_flags,
                        int m, const double *A, double *result )
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dplssq(
			quark, task_flags,
            m, A, result );
	} 
	else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task inout([m]A) out([1]result) label(dplssq)
		CORE_dplssq(m, A, result);
	}
}
