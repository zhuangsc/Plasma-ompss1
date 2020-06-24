/**
 *
 * @file pdlaset.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, double, m, n)
/***************************************************************************//**
 **/
void plasma_pdsrtdg_quark( PLASMA_desc A, double *work, 
                              PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int i, j;
    int ldai, ldaj;
    int tempim;
    int tempjm, tempjn;
    int minmn = min(A.mt, A.nt);

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;

    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

       for (j = 0; j < A.mt; j++){
           tempjm = j == A.mt-1 ? A.m-j*A.mb : A.mb;
           tempjn = j == A.nt-1 ? A.n-j*A.nb : A.nb;
           ldaj = BLKLDD(A, j);
           RT_dsrtdg( plasma->quark, &task_flags,
                 tempjm, A(j,j), ldaj, &work[j*A.nb], A.n);

       }
       RT_sort( A.n, work); 
}
