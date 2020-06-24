#include "runtime.h"

void RT_CORE_dswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, double *Aij,
                              int i1,  int i2, const int *ipiv, int inc,
                              const double *Akk, int ldak)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_dswptr_ontile(
			quark, task_flags,
			descA, Aij, i1, i2, ipiv, inc,
			Akk, ldak);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
                #pragma omp target device (smp) copy_deps
		#pragma omp task inout([1]Aij) in([(i2-i1+1)*abs(inc)]ipiv, [1]Akk) label(dswptr_ontile)
		CORE_dswptr_ontile(descA, i1, i2, ipiv, inc, Akk, ldak);
	}
}
