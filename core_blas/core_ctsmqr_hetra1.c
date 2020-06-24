/**
 *
 * @file core_ctsmqr_hetra1.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated c Tue Jan  7 11:44:48 2014
 *
 **/
#include <lapacke.h>
#include "common.h"
#undef REAL
#define COMPLEX

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_ctsmqr_hetra1: see CORE_ctsmqr
 *
 *  This kernel applies a left transformation on | A1'|
 *                                               | A2 |
 *
 * Needs therefore to make the explicit transpose of A1 before
 * and after the application of the block of reflectors
 * Can be further optimized by changing accordingly the underneath
 * kernel ztsrfb!
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg PlasmaLeft  : apply Q or Q**H from the Left;
 *         @arg PlasmaRight : apply Q or Q**H from the Right.
 *
 * @param[in] trans
 *         @arg PlasmaNoTrans   :  No transpose, apply Q;
 *         @arg PlasmaConjTrans :  ConjTranspose, apply Q**H.
 *
 * @param[in] m1
 *         The number of rows of the tile A1. M1 >= 0.
 *
 * @param[in] n1
 *         The number of columns of the tile A1. N1 >= 0.
 *
 * @param[in] m2
 *         The number of rows of the tile A2. M2 >= 0.
 *         M2 = M1 if side == PlasmaRight.
 *
 * @param[in] n2
 *         The number of columns of the tile A2. N2 >= 0.
 *         N2 = N1 if side == PlasmaLeft.
 *
 * @param[in] k
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] ib
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M1-by-N1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] lda1
 *         The leading dimension of the array A1. LDA1 >= max(1,M1).
 *
 * @param[in,out] A2
 *         On entry, the M2-by-N2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] lda2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M2).
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_CTSQRT in the first k columns of its array argument V.
 *
 * @param[in] ldv
 *         The leading dimension of the array V. LDV >= max(1,K).
 *
 * @param[in] T
 *         The IB-by-N1 triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] WORK
 *         Workspace array of size
 *             LDWORK-by-N1 if side == PlasmaLeft
 *             LDWORK-by-IB if side == PlasmaRight
 *
 * @param[in] ldwork
 *         The leading dimension of the array WORK.
 *             LDWORK >= max(1,IB) if side == PlasmaLeft
 *             LDWORK >= max(1,M1) if side == PlasmaRight
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ctsmqr_hetra1 = PCORE_ctsmqr_hetra1
#define CORE_ctsmqr_hetra1 PCORE_ctsmqr_hetra1
#define CORE_ctsmqr PCORE_ctsmqr
int  CORE_ctsmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M1, int N1, int M2, int N2, int K, int IB,
                 PLASMA_Complex32_t *A1, int LDA1,
                 PLASMA_Complex32_t *A2, int LDA2,
                 const PLASMA_Complex32_t *V, int LDV,
                 const PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *WORK, int LDWORK);
#endif
int CORE_ctsmqr_hetra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        PLASMA_Complex32_t *A1, int lda1,
                        PLASMA_Complex32_t *A2, int lda2,
                        const PLASMA_Complex32_t *V, int ldv,
                        const PLASMA_Complex32_t *T, int ldt,
                        PLASMA_Complex32_t *WORK, int ldwork)
{
    int i, j;

    if ( (m1 != n1) ) {
        coreblas_error(3, "Illegal value of M1, N1");
        return -3;
    }

    /* in-place transposition of A1 */
    for (j = 0; j < n1; j++){
        A1[j + j*lda1] = conjf(A1[j + j*lda1]);

        for (i = j+1; i < m1; i++){
            *WORK = *(A1 + i + j*lda1);
            *(A1 + i + j*lda1) = conjf(*(A1 + j + i*lda1));
            *(A1 + j + i*lda1) = conjf(*WORK);
        }
    }

    CORE_ctsmqr(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);

    /* in-place transposition of A1 */
    for (j = 0; j < n1; j++){
        A1[j + j*lda1] = conjf(A1[j + j*lda1]);

        for (i = j+1; i < m1; i++){
            *WORK = *(A1 + i + j*lda1);
            *(A1 + i + j*lda1) = conjf(*(A1 + j + i*lda1));
            *(A1 + j + i*lda1) = conjf(*WORK);
        }
    }

    return PLASMA_SUCCESS;
}
