#include "runtime.h"

#if 0
/* Maximum number of simultaneous tasks */
#define MAX_TASK 256
struct rt_s_tasks{
	int n_task;
	volatile int *task_pool;
} RT_TASK;
#endif

int RT_get_runtime()
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	return plasma->runtime;
}

void RT_get_size()
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if ( RT_DEBUG ) {
		fprintf(OUTFILE, "PLASMA_SIZE: %d\n", PLASMA_SIZE);
	}
}

int RT_get_ws()
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	return PLASMA_SIZE;
}

void RT_set_ws(int sze)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	PLASMA_SIZE = sze;
//	if ( RT_DEBUG ) {
//		fprintf(OUTFILE, "change PLASMA_SIZE to %d\n", sze);
//	}
}

void RT_runtime_info()
{
	if ( RT_DEBUG ) {
		plasma_context_t *plasma;
		plasma = plasma_context_self();
		fprintf(OUTFILE, "PLASMA_SIZE: %d ", PLASMA_SIZE);
		switch (plasma->runtime) 
		{
			case PLASMA_QUARK:
				fprintf(OUTFILE, "QUARK in use ");
				break;
			case PLASMA_OMPSS:
				fprintf(OUTFILE, "OmpSs in use ");
				break;
			default:
				fprintf(OUTFILE, "Unknown scheduler ");
		}
		switch (plasma->scheduling)
		{
			case PLASMA_STATIC_SCHEDULING:
				fprintf(OUTFILE, "Static scheduling ");
				break;
			case PLASMA_DYNAMIC_SCHEDULING:
				fprintf(OUTFILE, "Dynamic scheduling ");
				break;
			default:
				fprintf(OUTFILE, " ");
		}
		fprintf(OUTFILE, "\n");
	}
}

void RT_dynamic_sync()
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();

	if (plasma->runtime == PLASMA_QUARK) {
		plasma_dynamic_sync();
	}
	else if(plasma->runtime == PLASMA_OMPSS) {
		//#pragma omp taskwait noflush
		#pragma omp taskwait
	}
}

void RT_dynamic_sync_on(double *ptr)
{
	#pragma omp taskwait on(ptr)
}

void RT_CORE_free(Quark *quark, Quark_Task_Flags *task_flags, void *A, int szeA)
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	if (plasma->runtime == PLASMA_QUARK) {
		QUARK_CORE_free(quark, task_flags, A, szeA);
	}
	else if (plasma->runtime == PLASMA_OMPSS) {
		if (A != NULL) {
			printf("");
//			#pragma omp task inout([szeA]A)
//			free(A);
		}
	}
}

#pragma omp task inout([sze]A)
void RT_CORE_foo(double *A, int sze)
{
	return;
}

#pragma omp task concurrent([szeA]A) inout([szeB]B)
void RT_CORE_foo2(double *A, int szeA, double *B, int szeB)
{
	return;
}

int RT_threads_num()
{
	int threads_num = omp_get_num_threads();
	return threads_num;
}

int RT_global_rank()
{
	int grank = omp_get_thread_num();
	return grank;
}

#if 0
void RT_init()
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();

	if ( plasma->runtime == PLASMA_OMPSS ) {
		RT_TASK.n_task = 0;
		RT_TASK.task_pool = malloc(MAX_TASK * sizeof(int));
		memset(RT_TASK.task_pool, -1, sizeof(int)*MAX_TASK);
	}
}

int RT_task_init()
{
	plasma_context_t *plasma;
	plasma = plasma_context_self();
	int id = 0;
	if ( plasma->runtime == PLASMA_OMPSS ) {
		int *pool = RT_TASK.task_pool;
		#pragma omp critical 
		{
			while ( pool[id] != -1 ) {
				id++;
			}
			if ( id <= MAX_TASK ) {
				RT_TASK.n_task++;
				pool[id] = 0;
			} else {
				fprintf(stderr, "RT_TASK full\n");
			}
		}
	}
	return id;
}

int RT_local_rank(int t_id)
{
	int grank = RT_global_rank();
	int lrank;
	#pragma omp critical
	{
		lrank = RT_TASK.task_pool[t_id];
		RT_TASK.task_pool[t_id]++;
		printf("t_id: %d, grank: %d, lrank: %d\n", t_id, grank, lrank);
	}
	return lrank;
}

void RT_task_fini(int t_id)
{
	#pragma omp critical
	{
		RT_TASK.task_pool[t_id] = -1;
		RT_TASK.n_task--;
	}
}

void RT_fini()
{
	RT_TASK.n_task = 0;
	free(RT_TASK.task_pool);
}
#endif
