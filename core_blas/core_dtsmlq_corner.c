/**
 *
 * @file core_dtsmlq_corner.c
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
 *  CORE_dtsmlq_corner: see CORE_dtsmlq
 *
 * This kernel applies left and right transformations as depicted below:
 * |I -VTV'| * | A1  A2 | * |I - VT'V'|
 *             | A2' A3 |
 * where A1 and A3 are symmetric matrices.
 * Only the lower part is referenced.
 * This is an adhoc implementation, can be further optimized...
 *
 *******************************************************************************
 *
 * @param[in] m1
 *         The number of rows of the tile A1. m1 >= 0.
 *
 * @param[in] n1
 *         The number of columns of the tile A1. n1 >= 0.
 *
 * @param[in] m2
 *         The number of rows of the tile A2. m2 >= 0.
 *
 * @param[in] n2
 *         The number of columns of the tile A2. n2 >= 0.
 *
 * @param[in] m3
 *         The number of rows of the tile A3. m3 >= 0.
 *
 * @param[in] n3
 *         The number of columns of the tile A3. n3 >= 0.
 *
 * @param[in] k
 *         The number of elementary reflectors whose product defines
 *         the matrix Q.
 *
 * @param[in] ib
 *         The inner-blocking size.  ib >= 0.
 *
 * @param[in] nb
 *         The blocking size.  nb >= 0.
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
 * @param[in,out] A3
 *         On entry, the m3-by-n3 tile A3.
 *
 * @param[in] lda3
 *         The leading dimension of the tile A3. lda3 >= max(1,m3).
 *
 * @param[in] V
 *         The i-th row must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_DTSLQT in the first k rows of its array argument V.
 *
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
#pragma weak CORE_dtsmlq_corner = PCORE_dtsmlq_corner
#define CORE_dtsmlq_corner PCORE_dtsmlq_corner
#define CORE_dtsmlq PCORE_dtsmlq
int  CORE_dtsmlq(PLASMA_enum side, PLASMA_enum trans,
                 int m1, int n1, int m2, int n2, int K, int IB,
                 double *A1, int lda1,
                 double *A2, int lda2,
                 const double *V, int ldv,
                 const double *T, int ldt,
                 double *WORK, int LDWORK);
#endif
int CORE_dtsmlq_corner( int m1, int n1, int m2, int n2, int m3, int n3,
                        int k, int ib, int nb,
                        double *A1, int lda1,
                        double *A2, int lda2,
                        double *A3, int lda3,
                        const double *V, int ldv,
                        const double *T, int ldt,
                        double *WORK, int ldwork)
{
    PLASMA_enum side;
    PLASMA_enum trans;
    int i, j;

    if ( m1 != n1 ) {
        coreblas_error(1, "Illegal value of M1, N1");
        return -1;
    }

    /* Rebuild the symmetric block: WORK <- A1 */
    for (i = 0; i < m1; i++)
        for (j = i; j < n1; j++){
            *(WORK + i + j*ldwork) = *(A1 + i + j*lda1);
            if (j > i){
                *(WORK + j + i*ldwork) =  ( *(WORK + i + j*ldwork) );
            }
        }

    /*  Copy the transpose of A2: WORK+nb*ldwork <- A2' */
    for (j = 0; j < n2; j++)
        for (i = 0; i < m2; i++){
            *(WORK + j + (i + nb) * ldwork) = ( *(A2 + i + j*lda2) );
        }

    side = PlasmaRight;
    trans = PlasmaTrans;

    /*  Right application on |A1 A2| */
    CORE_dtsmlq(side, trans, m1, n1, m2, n2, k, ib,
                WORK, ldwork, A2, lda2,
                V, ldv, T, ldt,
                WORK+3*nb*ldwork, ldwork);

    /*  Rebuild the symmetric block: WORK+2*nb*ldwork <- A3 */
    for (i = 0; i < m3; i++)
        for (j = i; j < n3; j++){
            *(WORK + i + (j + 2*nb) * ldwork) = *(A3 + i + j*lda3);
            if (j > i){
                *(WORK + j + (i + 2*nb) * ldwork) =   ( *(WORK + i + (j + 2*nb) * ldwork) );
            }
        }

    /*  Right application on | A2' A3 | */
    CORE_dtsmlq(side, trans, n2, m2, m3, n3, k, ib,
                WORK+nb*ldwork, ldwork, WORK+2*nb*ldwork, ldwork,
                V, ldv, T, ldt,
                WORK + 3*nb*ldwork, ldwork);

    side = PlasmaLeft;
    trans = PlasmaNoTrans;

    /*  Left application on | A1  | */
    /*                      | A2' | */
    CORE_dtsmlq(side, trans, m1, n1, n2, m2, k, ib,
                WORK, ldwork, WORK+nb*ldwork, ldwork,
                V, ldv, T, ldt,
                WORK + 3*nb*ldwork, ldwork);

    /*  Copy back the final result to the upper part of A1 */
    /*  A1 = WORK */
    for (i = 0; i < m1; i++)
        for (j = i; j < n1; j++)
            *(A1 + i + j*lda1) = *(WORK + i + j*ldwork);

    /*  Left application on | A2 | */
    /*                      | A3 | */
    CORE_dtsmlq(side, trans, m2, n2, m3, n3, k, ib,
                A2, lda2, WORK+2*nb*ldwork, ldwork,
                V, ldv, T, ldt,
                WORK + 3*nb*ldwork, ldwork);

    /*  Copy back the final result to the upper part of A3 */
    /*  A3 = WORK+2*nb*ldwork */
    for (i = 0; i < m3; i++)
        for (j = i; j < n3; j++)
            *(A3 + i + j*lda3) = *(WORK + i + (j+ 2*nb) * ldwork);

    return PLASMA_SUCCESS;
}
