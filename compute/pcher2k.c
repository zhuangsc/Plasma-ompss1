/**
 *
 * @file pcher2k.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Tue Jan  7 11:45:12 2014
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex32_t, m, n)
#define B(m,n) BLKADDR(B, PLASMA_Complex32_t, m, n)
#define C(m,n) BLKADDR(C, PLASMA_Complex32_t, m, n)
/***************************************************************************//**
 *  Parallel tile Hermitian rank-k update - static scheduling
 **/
void plasma_pcher2k(plasma_context_t *plasma)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    PLASMA_Complex32_t alpha;
    PLASMA_desc A;
    PLASMA_desc B;
    float beta;
    PLASMA_desc C;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int m, n, k;
    int next_m;
    int next_n;
    int ldam, ldan, ldak;
    int ldbm, ldbn, ldbk;
    int ldcm, ldcn;
    int tempkn, tempkm, tempmm, tempnn;

    PLASMA_Complex32_t zone = (PLASMA_Complex32_t)1.0;
    PLASMA_Complex32_t zbeta;
    float dbeta;

    plasma_unpack_args_9(uplo, trans, alpha, A, B, beta, C, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= C.mt && n < C.nt) {
        n++;
        m = m-C.mt+n;
    }

    while (n < C.nt) {
        next_n = n;
        next_m = m + PLASMA_SIZE;
        while (next_m >= C.mt && next_n < C.nt) {
            next_n++;
            next_m = next_m - C.mt + next_n;
        }

        tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
        tempnn = n == C.nt-1 ? C.n-n*C.nb : C.nb;

        ldcn = BLKLDD(C, n);
        ldcm = BLKLDD(C, m);

        if (m == n) {
            /*
             *  PlasmaNoTrans
             */
            if (trans == PlasmaNoTrans) {
                ldam = BLKLDD(A, m);
                ldbm = BLKLDD(B, m);
                for (k = 0; k < A.nt; k++) {
                    tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                    dbeta = k == 0 ? beta : 1.0;
                    CORE_cher2k(
                        uplo, trans,
                        tempnn, tempkn,
                        alpha, A(m, k), ldam,
                               B(m, k), ldbm,
                        dbeta, C(m, m), ldcm);
                }
            }
            /*
             *  Plasma[Conj]Trans
             */
            else {
                for (k = 0; k < A.mt; k++) {
                    tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    dbeta = k == 0 ? beta : 1.0;
                    CORE_cher2k(
                        uplo, trans,
                        tempnn, tempkm,
                        alpha, A(k, m), ldak,
                               B(k, m), ldbk,
                        dbeta, C(m, m), ldcm);
                }
            }
        }
        else {
            if (trans == PlasmaNoTrans) {
                ldam = BLKLDD(A, m);
                ldan = BLKLDD(A, n);
                ldbm = BLKLDD(B, m);
                ldbn = BLKLDD(B, n);
                /*
                 *  PlasmaNoTrans / PlasmaLower
                 */
                if (uplo == PlasmaLower) {
                    for (k = 0; k < A.nt; k++) {
                        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                        zbeta = k == 0 ? (PLASMA_Complex32_t)beta : zone;
                        CORE_cgemm(
                            trans, PlasmaConjTrans,
                            tempmm, tempnn, tempkn,
                            alpha, A(m, k), ldam,
                                   B(n, k), ldbn,
                            zbeta, C(m, n), ldcm);

                        CORE_cgemm(
                            trans, PlasmaConjTrans,
                            tempmm, tempnn, tempkn,
                            conjf(alpha), B(m, k), ldbm,
                                         A(n, k), ldan,
                            zone,        C(m, n), ldcm);
                    }
                }
                /*
                 *  PlasmaNoTrans / PlasmaUpper
                 */
                else {
                    for (k = 0; k < A.nt; k++) {
                        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                        zbeta = k == 0 ? (PLASMA_Complex32_t)beta : zone;
                        CORE_cgemm(
                            trans, PlasmaConjTrans,
                            tempnn, tempmm, tempkn,
                            alpha, A(n, k), ldan,
                                   B(m, k), ldbm,
                            zbeta, C(n, m), ldcn);

                        CORE_cgemm(
                            trans, PlasmaConjTrans,
                            tempnn, tempmm, tempkn,
                            conjf(alpha), B(n, k), ldbn,
                                         A(m, k), ldam,
                            zone,        C(n, m), ldcn);
                    }
                }
            }
            else {
                /*
                 *  Plasma[Conj]Trans / PlasmaLower
                 */
                if (uplo == PlasmaLower) {
                    for (k = 0; k < A.mt; k++) {
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                        zbeta = k == 0 ? (PLASMA_Complex32_t)beta : zone;
                        CORE_cgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm,
                            alpha, A(k, m), ldak,
                                   B(k, n), ldbk,
                            zbeta, C(m, n), ldcm);

                        CORE_cgemm(
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm,
                            conjf(alpha), B(k, m), ldbk,
                                         A(k, n), ldak,
                            zone,        C(m, n), ldcm);
                    }
                }
                /*
                 *  Plasma[Conj]Trans / PlasmaUpper
                 */
                else {
                    for (k = 0; k < A.mt; k++) {
                        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? (PLASMA_Complex32_t)beta : zone;
                        CORE_cgemm(
                            trans, PlasmaNoTrans,
                            tempnn, tempmm, tempkm,
                            alpha, A(k, n), ldak,
                                   B(k, m), ldbk,
                            zbeta, C(n, m), ldcm);

                        CORE_cgemm(
                            trans, PlasmaNoTrans,
                            tempnn, tempmm, tempkm,
                            conjf(alpha), B(k, n), ldbk,
                                         A(k, m), ldak,
                            zone,        C(n, m), ldcn);
                    }
                }
            }
        }
        m = next_m;
        n = next_n;
    }
}

/***************************************************************************//**
 *  Parallel tile Hermitian rank-k update - dynamic scheduling
 **/
void plasma_pcher2k_quark(PLASMA_enum uplo, PLASMA_enum trans,
                          PLASMA_Complex32_t alpha, PLASMA_desc A, PLASMA_desc B,
                          float beta,  PLASMA_desc C,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int m, n, k;
    int ldak, ldam, ldan, ldcm, ldcn;
    int ldbk, ldbm, ldbn;
    int tempnn, tempmm, tempkn, tempkm;

    PLASMA_Complex32_t zone   = (PLASMA_Complex32_t)1.0;
    PLASMA_Complex32_t zbeta;
    float dbeta;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (n = 0; n < C.nt; n++) {
        tempnn = n == C.nt-1 ? C.n-n*C.nb : C.nb;
        ldan = BLKLDD(A, n);
        ldbn = BLKLDD(B, n);
        ldcn = BLKLDD(C, n);
        /*
         *  PlasmaNoTrans
         */
        if (trans == PlasmaNoTrans) {
            for (k = 0; k < A.nt; k++) {
                tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                dbeta = k == 0 ? beta : 1.0;
                QUARK_CORE_cher2k(
                    plasma->quark, &task_flags,
                    uplo, trans,
                    tempnn, tempkn, A.mb,
                    alpha, A(n, k), ldan, /* ldan * K */
                           B(n, k), ldbn,
                    dbeta, C(n, n), ldcn); /* ldc  * N */
            }
            /*
             *  PlasmaNoTrans / PlasmaLower
             */
            if (uplo == PlasmaLower) {
                for (m = n+1; m < C.mt; m++) {
                    tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    ldcm = BLKLDD(C, m);
                    for (k = 0; k < A.nt; k++) {
                        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                        zbeta = k == 0 ? (PLASMA_Complex32_t)beta : zone;
                        QUARK_CORE_cgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaConjTrans,
                            tempmm, tempnn, tempkn, A.mb,
                            alpha, A(m, k), ldam,  /* ldam * K */
                                   B(n, k), ldbn,  /* ldan * K */
                            zbeta, C(m, n), ldcm); /* ldc  * N */

                        QUARK_CORE_cgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaConjTrans,
                            tempmm, tempnn, tempkn, A.mb,
                            conjf(alpha), B(m, k), ldbm,  /* ldam * K */
                                   A(n, k), ldan,  /* ldan * K */
                            zone,  C(m, n), ldcm); /* ldc  * N */
                    }
                }
            }
            /*
             *  PlasmaNoTrans / PlasmaUpper
             */
            else {
                for (m = n+1; m < C.mt; m++) {
                    tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    for (k = 0; k < A.nt; k++) {
                        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                        zbeta = k == 0 ? (PLASMA_Complex32_t)beta : zone;
                        QUARK_CORE_cgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaConjTrans,
                            tempnn, tempmm, tempkn, A.mb,
                            alpha, A(n, k), ldan,  /* ldan * K */
                                   B(m, k), ldbm,  /* ldam * M */
                            zbeta, C(n, m), ldcn); /* ldc  * M */

                        QUARK_CORE_cgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaConjTrans,
                            tempnn, tempmm, tempkn, A.mb,
                            conjf(alpha), B(n, k), ldan,  /* ldan * K */
                                   A(m, k), ldam,  /* ldam * M */
                            zone,  C(n, m), ldcn); /* ldc  * M */
                    }
                }
            }
        }
        /*
         *  Plasma[Conj]Trans
         */
        else {
            for (k = 0; k < A.mt; k++) {
                tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);
                dbeta = k == 0 ? beta : 1.0;
                QUARK_CORE_cher2k(
                    plasma->quark, &task_flags,
                    uplo, trans,
                    tempnn, tempkm, A.mb,
                    alpha, A(k, n), ldak,  /* lda * N */
                           B(k, n), ldbk,
                    dbeta, C(n, n), ldcn); /* ldc * N */
            }
            /*
             *  Plasma[Conj]Trans / PlasmaLower
             */
            if (uplo == PlasmaLower) {
                for (m = n+1; m < C.mt; m++) {
                    tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
                    ldcm = BLKLDD(C, m);
                    for (k = 0; k < A.mt; k++) {
                        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? (PLASMA_Complex32_t)beta : zone;
                        QUARK_CORE_cgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm, A.mb,
                            alpha, A(k, m), ldak,  /* lda * M */
                                   B(k, n), ldbk,  /* lda * N */
                            zbeta, C(m, n), ldcm); /* ldc * N */

                        QUARK_CORE_cgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaNoTrans,
                            tempmm, tempnn, tempkm, A.mb,
                            conjf(alpha), B(k, m), ldbk,  /* lda * M */
                                   A(k, n), ldak,  /* lda * N */
                            zone,  C(m, n), ldcm); /* ldc * N */
                    }
                }
            }
            /*
             *  Plasma[Conj]Trans / PlasmaUpper
             */
            else {
                for (m = n+1; m < C.mt; m++) {
                    tempmm = m == C.mt-1 ? C.m-m*C.mb : C.mb;
                    for (k = 0; k < A.mt; k++) {
                        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                        ldak = BLKLDD(A, k);
                        ldbk = BLKLDD(B, k);
                        zbeta = k == 0 ? (PLASMA_Complex32_t)beta : zone;
                        QUARK_CORE_cgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaNoTrans,
                            tempnn, tempmm, tempkm, A.mb,
                            alpha, A(k, n), ldak,  /* lda * K */
                                   B(k, m), ldbk,  /* lda * M */
                            zbeta, C(n, m), ldcn); /* ldc * M */

                        QUARK_CORE_cgemm(
                            plasma->quark, &task_flags,
                            trans, PlasmaNoTrans,
                            tempnn, tempmm, tempkm, A.mb,
                            conjf(alpha), B(k, n), ldbk,  /* lda * K */
                                   A(k, m), ldak,  /* lda * M */
                            zone,  C(n, m), ldcn); /* ldc * M */
                    }
                }
            }
        }
    }
}
