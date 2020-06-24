/**
 *
 * @file psplgsy.c
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
#include "common.h"

#define A(m,n) BLKADDR(A, float, m, n)
/***************************************************************************//**
 *  Parallel tile Cholesky factorization - static scheduling
 **/
void plasma_psplgsy(plasma_context_t *plasma)
{
    float bump;
    PLASMA_desc A;
    unsigned long long int seed;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int m, n;
    int next_m;
    int next_n;
    int ldam;
    int tempmm, tempnn;

    plasma_unpack_args_5(bump, A, seed, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= A.mt) {
        n++;
        m = m - A.mt;
    }

    while ( n < A.nt ) {
        next_n = n;
        next_m = m;

        next_m += PLASMA_SIZE;
        while ( next_m >= A.mt && next_n < A.nt ) {
            next_n++;
            next_m = next_m - A.mt;
        }

        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
        ldam = BLKLDD(A, m);

        CORE_splgsy( 
            bump, tempmm, tempnn, A(m, n), ldam,
            A.m, m*A.mb, n*A.nb, seed );

        m = next_m;
        n = next_n;
    }
}

/***************************************************************************//**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void plasma_psplgsy_quark( float bump, PLASMA_desc A, unsigned long long int seed,
                           PLASMA_sequence *sequence, PLASMA_request *request )
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int m, n;
    int ldam;
    int tempmm, tempnn;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        ldam = BLKLDD(A, m);
        
        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            
            QUARK_CORE_splgsy( 
                plasma->quark, &task_flags,
                bump, tempmm, tempnn, A(m, n), ldam,
                A.m, m*A.mb, n*A.nb, seed );
        }
    }
}
