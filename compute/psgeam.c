/**
 *
 * @file psgeam.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 **/
#include "common.h"

#define A(m, n) BLKADDR(A, float, m, n)
#define B(m, n) BLKADDR(B, float, m, n)
void plasma_psgeam(plasma_context_t *plasma)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    float alpha;
    PLASMA_desc A;
    float beta;
    PLASMA_desc B;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int K, X, Y;
    int k, m, n;
    int next_m;
    int next_n;
    int ldam, ldbm, ldan, ldbn;
    float zbeta;
    float zone = (float)1.0;

    plasma_unpack_args_9(transA, transB, alpha, A, A, beta, B, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= A.mt && n < A.nt) {
        n++;
        m = m-A.mt;
    }
    while (n < A.nt) {
        next_m = m;
        next_n = n;

        next_m += PLASMA_SIZE;
        while (next_m >= A.mt && next_n < A.nt) {
            next_n++;
            next_m = next_m - A.mt;
        }

        X = m == A.mt-1 ? A.m - m*A.mb : A.mb;
        Y = n == A.nt-1 ? A.n - n*A.nb : A.nb;

        ldam = BLKLDD(A, m);
        ldbm = BLKLDD(B, m);
        /*
         *  A: PlasmaNoTrans / B: PlasmaNoTrans
         */
        if (transA == PlasmaNoTrans) {
            if (transB == PlasmaNoTrans) {
                    CORE_sgeam(
                        transA, transB,
                        X, Y,
                        alpha, A(m, n), ldam,
                        beta, B(m, n), ldbm);
            }
       /*
        *  A: PlasmaNoTrans / B: Plasma[Conj]Trans
        */
            else {
                ldbn = BLKLDD(B, n);
                    CORE_sgeam(
                        transA, transB,
                        X, Y,
                        alpha, A(m, n), ldam,
                        beta, B(n, m), ldbn);
            }
        }
        /*
         *  A: Plasma[Conj]Trans / B: PlasmaNoTrans
         */
        else {
            if (transB == PlasmaNoTrans) {
                    CORE_sgeam(
                        transA, transB,
                        X, Y,
                        alpha, A(n, m), ldan,
                        zbeta, B(m, n), ldbm);
            }
            /*
             *  A: Plasma[Conj]Trans / B: Plasma[Conj]Trans
             */
            else {
                ldbn = BLKLDD(B, n);
                    CORE_sgeam(
                        transA, transB,
                        X, Y,
                        alpha, A(n, m), ldan,
                        zbeta, B(n, m), ldbn);
            }
        }
        m = next_m;
        n = next_n;
    }
}


/***************************************************************************//**
 *  Parallel tile matrix-matrix multiplication - dynamic scheduling
 **/

void plasma_psgeam_quark(PLASMA_enum transA, PLASMA_enum transB,
                         float alpha, PLASMA_desc A, PLASMA_desc B,  
                         float beta, PLASMA_desc C,
                         PLASMA_sequence *sequence, PLASMA_request *request)
{   
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int m, n, k;
    int ldam, ldan, ldbn, ldbm;
    int tempmm, tempnn, tempkn, tempkm;

    float zbeta;
    float zone = (float)1.0;
    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    int X, Y;
    ldam = BLKLDD(A, m);
    ldbm = BLKLDD(B, m);
    for (m = 0; m < A.mt; m++) {
        X = m == A.mt-1 ? A.m-m*A.mb : A.mb;

        for (n = 0; n < A.nt; n++) {
            Y = n == A.nt-1 ? A.n-n*A.nb : A.nb;
             if (transA == PlasmaNoTrans) {
                if (transB == PlasmaNoTrans) {
                   QUARK_CORE_sgeam(
                      plasma->quark, &task_flags,
                      transA, transB,
                      X, Y, tempkm, A.mb,
                      alpha, A(m, n), ldam,
                      beta, B(m, n), ldbm);
                }
                else {
                   QUARK_CORE_sgeam(
                      plasma->quark, &task_flags,
                      transA, transB,
                      X, Y, tempkm, A.mb,
                      alpha, A(m, n), ldam,
                      beta, B(n, m), ldbm);
                }
              }
               else {
                    if (transB == PlasmaNoTrans) {
                       QUARK_CORE_sgeam(
                          plasma->quark, &task_flags,
                          transA, transB,
                          X, Y, tempkm, A.mb,
                          alpha, A(n, m), ldam,
                          beta, B(m, n), ldbm);
                     }
                     else {
                        QUARK_CORE_sgeam(
                        plasma->quark, &task_flags,
                        transA, transB,
                        X, Y, tempkm, A.mb,
                        alpha, A(n, m), ldam,
                        beta, B(n, m), ldbm);
                     }

             }

        }
    }
}
