/**
 *
 * @file psgeqrfrh.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:45:11 2014
 *
 **/
#include "common.h"

#define A(m,n)  BLKADDR(A, float, (m), (n))
#define T(m,n)  BLKADDR(T, float, (m), (n))
#define T2(m,n) BLKADDR(T, float, (m), (n)+A.nt)
/***************************************************************************//**
 *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
 **/
void plasma_psgeqrfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int M, RD;
    int ldaM, ldam, ldaMRD;
    int tempkmin, tempkn, tempMm, tempnn, tempmm, tempMRDm;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;

    for (k = 0; k < min(A.mt, A.nt); k++) {
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        for (M = k; M < A.mt; M += BS) {
            tempMm = M == A.mt-1 ? A.m-M*A.mb : A.mb;
            tempkmin = min(tempMm, tempkn);
            ldaM = BLKLDD(A, M);
            QUARK_CORE_sgeqrt(
                plasma->quark, &task_flags,
                tempMm, tempkn, ib, T.nb,
                A(M, k), ldaM,
                T(M, k), T.mb);
            for (n = k+1; n < A.nt; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_sormqr(
                    plasma->quark, &task_flags,
                    PlasmaLeft, PlasmaTrans,
                    tempMm, tempnn, tempkmin, ib, T.nb,
                    A(M, k), ldaM,
                    T(M, k), T.mb,
                    A(M, n), ldaM);
            }
            for (m = M+1; m < min(M+BS, A.mt); m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_stsqrt(
                    plasma->quark, &task_flags,
                    tempmm, tempkn, ib, T.nb,
                    A(M, k), ldaM,
                    A(m, k), ldam,
                    T(m, k), T.mb);

                for (n = k+1; n < A.nt; n++) {
                    tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                    QUARK_CORE_stsmqr(
                        plasma->quark, &task_flags,
                        PlasmaLeft, PlasmaTrans,
                        A.nb, tempnn, tempmm, tempnn, A.nb, ib, T.nb,
                        A(M, n), ldaM,
                        A(m, n), ldam,
                        A(m, k), ldam,
                        T(m, k), T.mb);
                }
            }
        }
        for (RD = BS; RD < A.mt-k; RD *= 2) {
            for (M = k; M+RD < A.mt; M += 2*RD) {
                tempMRDm = M+RD == A.mt-1 ? A.m-(M+RD)*A.mb : A.mb;
                ldaM   = BLKLDD(A, M   );
                ldaMRD = BLKLDD(A, M+RD);
                QUARK_CORE_sttqrt(
                    plasma->quark, &task_flags,
                    tempMRDm, tempkn, ib, T.nb,
                    A (M   , k), ldaM,
                    A (M+RD, k), ldaMRD,
                    T2(M+RD, k), T.mb);

                for (n = k+1; n < A.nt; n++) {
                    tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                    QUARK_CORE_sttmqr(
                        plasma->quark, &task_flags,
                        PlasmaLeft, PlasmaTrans,
                        A.nb, tempnn, tempMRDm, tempnn, A.nb, ib, T.nb,
                        A (M,    n), ldaM,
                        A (M+RD, n), ldaMRD,
                        A (M+RD, k), ldaMRD,
                        T2(M+RD, k), T.mb);
                }
            }
        }
    }
}
