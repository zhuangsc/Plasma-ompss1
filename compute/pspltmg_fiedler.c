/**
 *
 * @file pspltmg_fiedler.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:45:12 2014
 *
 **/
#include <stdlib.h>
#include "common.h"

#define A(m, n) BLKADDR(A, float, m, n)

/***************************************************************************//**
 *  Parallel tile Fiedler matrix generation - dynamic scheduling
 **/
void plasma_pspltmg_fiedler_quark( PLASMA_desc A, unsigned long long int seed,
                                   PLASMA_sequence *sequence, PLASMA_request *request )
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    float **work;
    int m, n;
    int ldam;
    int tempm0;
    int tempmm, tempnn;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    work = (float**) malloc ( A.mt * sizeof( float* ) );

    /* Generate a random vector of size A.m */
    for (m = 0; m < A.mt; m++) {
        tempm0 = m * A.mb;
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;

        /* Allocate temporary vector and initialize it randomly */
        work[m] = (float*)plasma_shared_alloc(plasma, tempmm, PlasmaRealFloat);

        QUARK_CORE_splrnt(
            plasma->quark, &task_flags,
            tempmm, 1, work[m], tempmm,
            A.m, tempm0+1, 0, seed );
    }

    /* Compute the fiedler Matrix tile by tile */
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        ldam = BLKLDD(A, m);

        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;

            QUARK_CORE_spltmg_fiedler(
                plasma->quark, &task_flags,
                tempmm, tempnn,
                work[m], 1,
                work[n], 1,
                A(m, n), ldam );
        }
    }

    /* Submit the workspace free */
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        QUARK_CORE_free(plasma->quark, &task_flags,
                        work[m], (tempmm)*sizeof(float));
    }

    /* We can free work and loose all pointers because they are already saved by Quark */
    free(work);
}
