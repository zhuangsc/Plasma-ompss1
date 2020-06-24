/**
 *
 * @file core_dtsmlq_sytra1.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:44:48 2014
 *
 **/
#include <lapacke.h>
#include "common.h"
#undef COMPLEX
#define REAL

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 * CORE_dtsmlq_sytra1: see CORE_dtsmlq
 *
 * This kernel applies a Right transformation on | A1' A2 |
 * and does not handle the transpose of A1.
 * Needs therefore to make the explicit transpose of A1 before
 * and after the application of the block of reflectors
 * Can be further optimized by changing accordingly the underneath
 * kernel ztsrfb!
 *
 *******************************************************************************
 *
 * @param[in] side
 *         @arg PlasmaLeft  : apply Q or Q**T from the Left;
 *         @arg PlasmaRight : apply Q or Q**T from the Right.
 *
 * @param[in] trans
 *         @arg PlasmaNoTrans   :  No transpose, apply Q;
 *         @arg PlasmaTrans :  ConjTranspose, apply Q**T.
 *
 * @param[in] m1
 *         The number of rows of the tile A1. m1 >= 0.
 *
 * @param[in] n1
 *         The number of columns of the tile A1. n1 >= 0.
 *
 * @param[in] m2
 *         The number of rows of the tile A2. m2 >= 0.
 *         m2 = m1 if side == PlasmaRight.
 *
 * @param[in] n2
 *         The number of columns of the tile A2. n2 >= 0.
 *         n2 = n1 if side == PlasmaLeft.
 *
 * @param[in] k
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] ib
 *         The inner-blocking size.  ib >= 0.
 *
 * @param[in,out] A1
 *         On entry, the m1-by-n1 tile A1.
 *         On exit, A1 is overwritten by the application of Q.
 *
 * @param[in] lda1
 *         The leading dimension of the array A1. lda1 >= max(1,m1).
 *
 * @param[in,out] A2
 *         On entry, the m2-by-n2 tile A2.
 *         On exit, A2 is overwritten by the application of Q.
 *
 * @param[in] lda2
 *         The leading dimension of the tile A2. lda2 >= max(1,m2).
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_DTSLQT in the first k rows of its array argument V.
 *!
 * @param[in] ldv
 *         The leading dimension of the array V. ldv >= max(1,K).
 *
 * @param[in] T
 *         The IB-by-n1 triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. ldt >= IB.
 *
 * @param[out] WORK
 *         Workspace array of size
 *             LDWORK-by-m1 if side == PlasmaLeft
 *             LDWORK-by-IB if side == PlasmaRight
 *
 * @param[in] ldwork
 *         The leading dimension of the array WORK.
 *             LDWORK >= max(1,IB) if side == PlasmaLeft
 *             LDWORK >= max(1,n1) if side == PlasmaRight
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtsmlq_sytra1 = PCORE_dtsmlq_sytra1
#define CORE_dtsmlq_sytra1 PCORE_dtsmlq_sytra1
#define CORE_dtsmlq PCORE_dtsmlq
int  CORE_dtsmlq(PLASMA_enum side, PLASMA_enum trans,
                 int m1, int n1, int m2, int n2, int K, int IB,
                 double *A1, int lda1,
                 double *A2, int lda2,
                 const double *V, int ldv,
                 const double *T, int ldt,
                 double *WORK, int LDWORK);
#endif
int CORE_dtsmlq_sytra1( PLASMA_enum side, PLASMA_enum trans,
                        int m1, int n1, int m2, int n2,
                        int k, int ib,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        const double *V, int ldv,
                        const double *T, int ldt,
                        double *WORK, int ldwork)
{
    int i, j;

    if ( (m1 != n1) ) {
        coreblas_error(3, "Illegal value of M1, N1");
        return -3;
    }

    /* in-place transposition of A1 */
    for (j = 0; j < n1; j++){
        A1[j + j*lda1] = (A1[j + j*lda1]);

        for (i = j+1; i < m1; i++){
            *WORK = *(A1 + i + j*lda1);
            *(A1 + i + j*lda1) = (*(A1 + j + i*lda1));
            *(A1 + j + i*lda1) = (*WORK);
        }
    }

    CORE_dtsmlq(side, trans, m1, n1, m2, n2, k, ib,
                A1, lda1, A2, lda2,
                V,  ldv,  T,  ldt,
                WORK, ldwork);

    /* in-place transposition of A1 */
    for (j = 0; j < n1; j++){
        A1[j + j*lda1] = (A1[j + j*lda1]);

        for (i = j+1; i < m1; i++){
            *WORK = *(A1 + i + j*lda1);
            *(A1 + i + j*lda1) = (*(A1 + j + i*lda1));
            *(A1 + j + i*lda1) = (*WORK);
        }
    }

    return PLASMA_SUCCESS;
}
