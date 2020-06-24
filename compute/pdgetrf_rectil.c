/**
 *
 * @file pdgetrf_rectil.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * LU with Partial pivoting.
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @author Hatem Ltaief
 * @date 2009-11-15
 *
 * @generated d Tue Jan  7 11:45:13 2014
 *
 **/
#include "common.h"

void CORE_dgetrf_rectil_init(void);

#define PARALLEL_KERNEL
#define A(m,n) BLKADDR(A, double, m, n)
#define IPIV(k) &(IPIV[(int64_t)A.mb*(int64_t)(k)])

#define plasma_pdgetrf_rectil_rl_quark plasma_pdgetrf_rectil_quark

/***************************************************************************//**
 *  Parallel tile LU factorization - dynamic scheduling - Right looking
 **/
void plasma_pdgetrf_rectil_rl_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request)
{
    int k, m, n;
    plasma_context_t *plasma;
    int tempk, tempm, tempkm, tempkn, tempmm, tempnn;
    int ldak, ldam;
    Quark_Task_Flags task_flagsP = Quark_Task_Flags_Initializer;
    Quark_Task_Flags task_flagsU = Quark_Task_Flags_Initializer;

    double zone  = (double)1.0;
    double mzone = (double)-1.0;

	/* OmpSs doesn't accept fake address as task dependencies*/
    double * fakedep0 = malloc(min(A.mt, A.nt) * sizeof(double));
	double * fakedep;

     /* How many threads per panel? Probably needs to be adjusted during factorization. */
    int panel_thread_count;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flagsP, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    QUARK_Task_Flag_Set(&task_flagsU, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    /* We start at PLASMA_SIZE-1, to keep the first thread adding task to the queue */
    panel_thread_count = min( max(PLASMA_SIZE, 2)-1, 48); /* kernel doesn't accept more than 48 cores */

    QUARK_Task_Flag_Set(&task_flagsP, TASK_THREAD_COUNT, panel_thread_count );

    CORE_dgetrf_rectil_init();

    for (k = 0; k < min(A.mt, A.nt); k++)
    {
        tempk  = k * A.mb;
        tempm  = A.m - tempk;
        tempkm = k == A.mt-1 ? tempm      : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        ldak = BLKLDD(A, k);

        QUARK_Task_Flag_Set(&task_flagsU, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - k );

        while ( ((panel_thread_count * 4 * A.mb) > tempm)
                && (panel_thread_count > 1) ) {
          panel_thread_count--;
          QUARK_Task_Flag_Set(&task_flagsP, TASK_THREAD_COUNT, panel_thread_count );
        }

        QUARK_Task_Flag_Set(&task_flagsP, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - k );

        RT_CORE_dgetrf_rectil(
            plasma->quark, &task_flagsP,
            plasma_desc_submatrix(A, tempk, k*A.nb, tempm, tempkn),
            A(k, k), A.mb*A.nb, IPIV(k),
            sequence, request, 1, tempk,
            panel_thread_count );

        /*
         * Update the trailing submatrix
         */
		fakedep = fakedep0 + k;
//        fakedep = (void *)(intptr_t)(k+1);
        for (n = k+1; n < A.nt; n++)
        {
            QUARK_Task_Flag_Set(&task_flagsU, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - n );
            /*
             * Apply row interchange after the panel (work on the panel)
             */
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            RT_CORE_dswptr_ontile(
                plasma->quark, &task_flagsU,
                plasma_desc_submatrix(A, tempk, n*A.nb, tempm, tempnn),
                A(k, n), 1, tempkm, IPIV(k), 1,
                A(k, k), ldak);

            m = k+1;
            if ( m < A.mt ) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);

                RT_CORE_dgemm2(
                    plasma->quark, &task_flagsU,
                    PlasmaNoTrans, PlasmaNoTrans,
                    tempmm, tempnn, A.nb, A.mb,
                    mzone, A(m, k), ldam,
                           A(k, n), ldak,
                    zone,  A(m, n), ldam);

                for (m = k+2; m < A.mt; m++)
                {
                    tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                    ldam = BLKLDD(A, m);

                    RT_CORE_dgemm_f2(
                        plasma->quark, &task_flagsU,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempmm, tempnn, A.nb, A.mb,
                        mzone, A(m, k), ldam,
                               A(k, n), ldak,
                        zone,  A(m, n), ldam,
                        /* Dependency on next swapa (gemm need to be done before) */
                        A(k+1, n), A.mb*A.nb, INOUT | GATHERV,
                        /* Dependency on next swapb (gemm need to use panel k before it has to be swaped */
                        fakedep,   1,         INPUT );
                }
            }
        }
    }

    QUARK_Task_Flag_Set(&task_flagsU, TASK_PRIORITY, QUARK_TASK_MIN_PRIORITY );
    for (k = 0; k < min(A.mt, A.nt); k++)
    {
        int mintmp;
        tempk  = k * A.mb;
        tempm  = A.m - tempk;
        tempkm = k == A.mt-1 ? tempm : A.mb;
        tempkn = k == A.nt-1 ? A.n - k * A.nb : A.nb;
        mintmp = min(tempkm, tempkn);
        ldak = BLKLDD(A, k);

        /*
         * Apply row interchange behind the panel (work on the panel)
         */
		fakedep = fakedep0 + k;
//        fakedep = (void*)(intptr_t)k;
        for (n = 0; n < k; n++)
        {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            RT_CORE_dlaswp_ontile_f2(
                plasma->quark, &task_flagsU,
                plasma_desc_submatrix(A, tempk, n*A.nb, tempm, tempnn),
                A(k, n), 1, mintmp, IPIV(k), 1,
                /* Dependency on previous swapb */
                A(k-1,n), A.lm*A.nb, INPUT,
                /* Dependency on all GEMM from previous step */
                fakedep,  1,         INOUT | GATHERV );
        }
    }
	free(fakedep0);
}

/***************************************************************************//**
 *  Parallel tile LU factorization - dynamic scheduling - Left looking
 **/
void plasma_pdgetrf_tile_ll_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request)
{
    int k, m, n;
    plasma_context_t *plasma;
    int tempkm, tempkn, tempmm, tempnn;
    int tempk, tempm;
    int ldak, ldam;
    Quark_Task_Flags task_flagsP = Quark_Task_Flags_Initializer;
    Quark_Task_Flags task_flagsU = Quark_Task_Flags_Initializer;

    double zone  = (double)1.0;
    double mzone = (double)-1.0;

	/* OmpSs doesn't accept fake address as task dependencies*/
    double * fakedep0 = malloc(min(A.mt, A.nt) * sizeof(double));
	double * fakedep;

     /* How many threads per panel? Probably needs to be adjusted during factorization. */
    int panel_thread_count;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flagsP, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    QUARK_Task_Flag_Set(&task_flagsU, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    /* We start at PLASMA_SIZE-1, to keep the first thread adding task to the queue */
    panel_thread_count = min( max(PLASMA_SIZE, 2)-1, 48); /* kernel doesn't accept more than 48 cores */

    QUARK_Task_Flag_Set(&task_flagsP, TASK_THREAD_COUNT, panel_thread_count );

    CORE_dgetrf_rectil_init();

	fakedep = fakedep0;
//    fakedep = (void*)(intptr_t)1;
    for (n = 0; n < A.nt; n++)
    {
        tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;

        QUARK_Task_Flag_Set(&task_flagsU, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - n );
        QUARK_Task_Flag_Set(&task_flagsP, TASK_PRIORITY, QUARK_TASK_MAX_PRIORITY - n );

        for (k = 0; k < min(A.mt, n); k++)
        {
            tempk  = k * A.mb;
            tempm  = A.m - tempk;
            tempkm = k == A.mt-1 ? tempm : A.mb;
            ldak = BLKLDD(A, k);

            RT_CORE_dswptr_ontile(
                plasma->quark, &task_flagsU,
                plasma_desc_submatrix(A, tempk, n*A.nb, tempm, tempnn),
                A(k, n), 1, tempkm, IPIV(k), 1,
                A(k, k), ldak);

            if (k < A.mt-1) {
                m = k+1;
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);

                RT_CORE_dgemm2(
                    plasma->quark, &task_flagsU,
                    PlasmaNoTrans, PlasmaNoTrans,
                    tempmm, tempnn, A.nb, A.mb,
                    mzone, A(m, k), ldam,
                           A(k, n), ldak,
                    zone,  A(m, n), ldam);

				fakedep = fakedep0 + k;
//                fakedep = (void*)(intptr_t)k;
                for (m = k+2; m < A.mt; m++)
                {
                    tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                    ldam = BLKLDD(A, m);

                    RT_CORE_dgemm_f2(
                        plasma->quark, &task_flagsU,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempmm, tempnn, A.nb, A.mb,
                        mzone, A(m, k), ldam,
                        A(k, n), ldak,
                        zone,  A(m, n), ldam,
                        /* Dependency on next swapa or getrf (gemm need to be done before) */
                        A(k+1, n), A.mb*A.nb, INOUT | GATHERV,
                        /* Dependency on next swapb (gemm need to use panel k before it has to be swaped */
                        fakedep,   1,         INPUT );
                }
            }
        }

        k = n;
        if ( n < A.mt ) {
            tempm  = A.m - k * A.mb;
            tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
            ldak = BLKLDD(A, k);

            while ( ((panel_thread_count * 4 * A.mb) > tempm)
                    && (panel_thread_count > 1) ) {
                panel_thread_count--;
                QUARK_Task_Flag_Set(&task_flagsP, TASK_THREAD_COUNT, panel_thread_count );
            }

            RT_CORE_dgetrf_rectil(
                plasma->quark, &task_flagsP,
                plasma_desc_submatrix(A, k*A.mb, k*A.nb, tempm, tempkn), A(k, k), A.mb*A.nb,
                IPIV(k), sequence, request,
                1, A.mb*k, panel_thread_count );
        }
    }

    QUARK_Task_Flag_Set(&task_flagsU, TASK_PRIORITY, QUARK_TASK_MIN_PRIORITY );
    for (k = 0; k < min(A.mt, A.nt); k++)
    {
        tempk  = k * A.mb;
        tempm  = A.m - tempk;
        tempkm = k == A.mt-1 ? tempm : A.mb;
        ldak = BLKLDD(A, k);

		fakedep = fakedep0 + k;
        fakedep = (void*)(intptr_t)k;
        for (n = 0; n < k; n++)
        {
            /*
             * Apply row interchange behind the panel (work on the panel)
             */
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            RT_CORE_dlaswp_ontile_f2(
                plasma->quark, &task_flagsU,
                plasma_desc_submatrix(A, tempk, n*A.nb, tempm, tempnn),
                A(k, n), 1, tempkm, IPIV(k), 1,
                /* Dependency on previous swapb */
                A(k-1,n), A.lm*A.nb, INPUT,
                /* Dependency on all GEMM from previous step */
                fakedep,  1,         INOUT | GATHERV );
        }
    }
	free(fakedep0);
}
