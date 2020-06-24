/**
 *
 * @file pstrtri.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:45:10 2014
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, float, m, n)
/***************************************************************************//**
 *  Parallel tile triangular matrix inverse - dynamic scheduling
 **/
void plasma_pstrtri_quark(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldam, ldan;
    int tempkn, tempmm, tempnn;

    float zone  = (float) 1.0;
    float mzone = (float)-1.0;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    /*
     *  PlasmaLower
     */
    if (uplo == PlasmaLower) {
        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            ldan = BLKLDD(A, n);
            for (m = n+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_strsm(
                    plasma->quark, &task_flags,
                    PlasmaRight, uplo, PlasmaNoTrans, diag,
                    tempmm, tempnn, A.mb,
                    mzone, A(n, n), ldan,
                           A(m, n), ldam);
            }
            for (m = n+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                for (k = 0; k < n; k++) {
                    tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                    QUARK_CORE_sgemm(
                        plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempmm, tempkn, tempnn, A.mb,
                        zone, A(m, n), ldam,
                              A(n, k), ldan,
                        zone, A(m, k), ldam);
                }
            }
            for (m = 0; m < n; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                QUARK_CORE_strsm(
                    plasma->quark, &task_flags,
                    PlasmaLeft, uplo, PlasmaNoTrans, diag,
                    tempnn, tempmm, A.mb,
                    zone, A(n, n), ldan,
                          A(n, m), ldan);
            }
            QUARK_CORE_strtri(
                plasma->quark, &task_flags,
                uplo, diag,
                tempnn, A.mb,
                A(n, n), ldan,
                sequence, request, A.nb*n);
        }
    }
    /*
     *  PlasmaUpper
     */
    else {
        for (m = 0; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            for (n = m+1; n < A.nt; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_strsm(
                    plasma->quark, &task_flags,
                    PlasmaLeft, uplo, PlasmaNoTrans, diag,
                    tempmm, tempnn, A.mb,
                    mzone, A(m, m), ldam,
                           A(m, n), ldam);
            }
            for (n = 0; n < m; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                ldan = BLKLDD(A, n);
                for (k = m+1; k < A.nt; k++) {
                    tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                    QUARK_CORE_sgemm(
                        plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempnn, tempkn, tempmm, A.mb,
                        zone, A(n, m), ldan,
                              A(m, k), ldam,
                        zone, A(n, k), ldan);
                }
                QUARK_CORE_strsm(
                    plasma->quark, &task_flags,
                    PlasmaRight, uplo, PlasmaNoTrans, diag,
                    tempnn, tempmm, A.mb,
                    zone, A(m, m), ldam,
                          A(n, m), ldan);
            }
            QUARK_CORE_strtri(
                plasma->quark, &task_flags,
                uplo, diag,
                tempmm, A.mb,
                A(m, m), ldam,
                sequence, request, A.mb*m);
        }
    }
}
