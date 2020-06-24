/**
 *
 * @file core_dpemv.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Dulceneia Becker
 * @date 2011-06-29
 * @generated d Tue Jan  7 11:44:49 2014
 *
 **/
#include <cblas.h>
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 *  CORE_dpemv  performs one of the matrix-vector operations
 *
 *     y = alpha*op( A )*x + beta*y
 *
 *  where  op( A ) is one of
 *
 *     op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**T,
 *
 *  alpha and beta are scalars, x and y are vectors and A is a
 *  pentagonal matrix (see further details).
 *
 *
 *  Arguments
 *  ==========
 *
 * @param[in] storev
 *
 *         @arg PlasmaColumnwise :  array A stored columwise
 *         @arg PlasmaRowwise    :  array A stored rowwise
 *
 * @param[in] trans
 *
 *         @arg PlasmaNoTrans   :  y := alpha*A*x    + beta*y.
 *         @arg PlasmaTrans     :  y := alpha*A**T*x + beta*y.
 *         @arg PlasmaTrans :  y := alpha*A**T*x + beta*y.
 *
 * @param[in] M
 *         Number of rows of the matrix A.
 *         M must be at least zero.
 *
 * @param[in] N
 *         Number of columns of the matrix A.
 *         N must be at least zero.
 *
 * @param[in] L
 *         Order of triangle within the matrix A (L specifies the shape
 *         of the matrix A; see further details).
 *
 * @param[in] ALPHA
 *         Scalar alpha.
 *
 * @param[in] A
 *         Array of size LDA-by-N.  On entry, the leading M by N part
 *         of the array A must contain the matrix of coefficients.
 *
 * @param[in] LDA
 *         Leading dimension of array A.
 *
 * @param[in] X
 *         On entry, the incremented array X must contain the vector x.
 *
 * @param[in] INCX
 *         Increment for the elements of X. INCX must not be zero.
 *
 * @param[in] BETA
 *         Scalar beta.
 *
 * @param[in,out] Y
 *         On entry, the incremented array Y must contain the vector y.
 *
 * @param[out] INCY
 *         Increment for the elements of Y. INCY must not be zero.
 *
 * @param[out] WORK
 *         Workspace array of size at least L.
 *
 *  Further Details
 *  ===============
 *
 *               |     N    |
 *            _   ___________   _
 *               |          |
 *     A:        |          |
 *          M-L  |          |
 *               |          |  M
 *            _  |.....     |
 *               \    :     |
 *            L    \  :     |
 *            _      \:_____|  _
 *
 *               |  L | N-L |
 *
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dpemv = PCORE_dpemv
#define CORE_dpemv PCORE_dpemv
#endif
int CORE_dpemv(PLASMA_enum trans, int storev,
               int M, int N, int L,
               double ALPHA,
               const double *A, int LDA,
               const double *X, int INCX,
               double BETA,
               double *Y, int INCY,
               double *WORK)
{

   /*
    *  y = alpha * op(A) * x + beta * y
    */

    int K;
    static double zzero = 0.0;


    /* Check input arguments */
    if ((trans != PlasmaNoTrans) && (trans != PlasmaTrans) && (trans != PlasmaTrans)) {
        coreblas_error(1, "Illegal value of trans");
        return -1;
    }
    if ((storev != PlasmaColumnwise) && (storev != PlasmaRowwise)) {
        coreblas_error(2, "Illegal value of storev");
        return -2;
    }
    if (!( ((storev == PlasmaColumnwise) && (trans != PlasmaNoTrans)) ||
           ((storev == PlasmaRowwise) && (trans == PlasmaNoTrans)) )) {
        coreblas_error(2, "Illegal values of trans/storev");
        return -2;
    }
    if (M < 0) {
        coreblas_error(3, "Illegal value of M");
        return -3;
    }
    if (N < 0) {
        coreblas_error(4, "Illegal value of N");
        return -4;
    }
    if (L > min(M ,N)) {
        coreblas_error(5, "Illegal value of L");
        return -5;
    }
    if (LDA < max(1,M)) {
        coreblas_error(8, "Illegal value of LDA");
        return -8;
    }
    if (INCX < 1) {
        coreblas_error(10, "Illegal value of INCX");
        return -10;
    }
    if (INCY < 1) {
        coreblas_error(13, "Illegal value of INCY");
        return -13;
    }

    /* Quick return */
    if ((M == 0) || (N == 0))
        return PLASMA_SUCCESS;
    if ((ALPHA == zzero) && (BETA == zzero))
        return PLASMA_SUCCESS;

    /* If L < 2, there is no triangular part */
    if (L == 1) L = 0;

    /* Columnwise */
    if (storev == PlasmaColumnwise) {
        /*
         *        ______________
         *        |      |     |    A1: A[ 0 ]
         *        |      |     |    A2: A[ M-L ]
         *        |  A1  |     |    A3: A[ (N-L) * LDA ]
         *        |      |     |
         *        |______| A3  |
         *        \      |     |
         *          \ A2 |     |
         *            \  |     |
         *              \|_____|
         *
         */


        /* Columnwise / NoTrans */
        if (trans == PlasmaNoTrans) {
            coreblas_error(1, "The case PlasmaNoTrans / PlasmaColumnwise is not yet implemented");
            return -1;
        }
        /* Columnwise / [Conj]Trans */
        else {

            /* L top rows of y */
            if (L > 0) {

                /* w = A_2' * x_2 */
                cblas_dcopy(
                            L, &X[INCX*(M-L)], INCX, WORK, 1);
                cblas_dtrmv(
                            CblasColMajor, (CBLAS_UPLO)PlasmaUpper,
                            (CBLAS_TRANSPOSE)trans,
                            (CBLAS_DIAG)PlasmaNonUnit,
                            L, &A[M-L], LDA, WORK, 1);

                if (M > L) {

                    /* y_1 = beta * y_1 [ + alpha * A_1 * x_1 ] */
                    cblas_dgemv(
                                CblasColMajor, (CBLAS_TRANSPOSE)trans,
                                M-L, L, (ALPHA), A, LDA,
                                X, INCX, (BETA), Y, INCY);

                    /* y_1 = y_1 + alpha * w */
                    cblas_daxpy(L, (ALPHA), WORK, 1, Y, INCY);

                } else {

                    /* y_1 = y_1 + alpha * w */
                    if (BETA == zzero) {
                        cblas_dscal(L, (ALPHA), WORK, 1);
                        cblas_dcopy(L, WORK, 1, Y, INCY);
                    } else {
                        cblas_dscal(L, (BETA), Y, INCY);
                        cblas_daxpy(L, (ALPHA), WORK, 1, Y, INCY);
                    }
                }
            }

            /* N-L bottom rows of Y */
            if (N > L) {
                K = N - L;
                cblas_dgemv(
                            CblasColMajor, (CBLAS_TRANSPOSE)trans,
                            M, K, (ALPHA), &A[LDA*L], LDA,
                            X, INCX, (BETA), &Y[INCY*L], INCY);
            }
        }
    }
    /* Rowwise */
    else {
        /*
         * --------------
         * |            | \           A1:  A[ 0 ]
         * |    A1      |   \         A2:  A[ (N-L) * LDA ]
         * |            | A2  \       A3:  A[ L ]
         * |--------------------\
         * |        A3          |
         * ----------------------
         *
         */

        /* Rowwise / NoTrans */
        if (trans == PlasmaNoTrans) {
            /* L top rows of A and y */
            if (L > 0) {

                /* w = A_2 * x_2 */
                cblas_dcopy(
                            L, &X[INCX*(N-L)], INCX, WORK, 1);
                cblas_dtrmv(
                            CblasColMajor, (CBLAS_UPLO)PlasmaLower,
                            (CBLAS_TRANSPOSE)PlasmaNoTrans,
                            (CBLAS_DIAG)PlasmaNonUnit,
                            L, &A[LDA*(N-L)], LDA, WORK, 1);

                if (N > L) {

                    /* y_1 = beta * y_1 [ + alpha * A_1 * x_1 ] */
                    cblas_dgemv(
                                CblasColMajor, (CBLAS_TRANSPOSE)PlasmaNoTrans,
                                L, N-L, (ALPHA), A, LDA,
                                X, INCX, (BETA), Y, INCY);

                    /* y_1 = y_1 + alpha * w */
                    cblas_daxpy(L, (ALPHA), WORK, 1, Y, INCY);

                } else {

                    /* y_1 = y_1 + alpha * w */
                    if (BETA == zzero) {
                        cblas_dscal(L, (ALPHA), WORK, 1);
                        cblas_dcopy(L, WORK, 1, Y, INCY);
                    } else {
                        cblas_dscal(L, (BETA), Y, INCY);
                        cblas_daxpy(L, (ALPHA), WORK, 1, Y, INCY);
                    }
                }
            }

            /* M-L bottom rows of Y */
            if (M > L) {
                cblas_dgemv(
                        CblasColMajor, (CBLAS_TRANSPOSE)PlasmaNoTrans,
                        M-L, N, (ALPHA), &A[L], LDA,
                        X, INCX, (BETA), &Y[INCY*L], INCY);
            }
        }
        /* Rowwise / [Conj]Trans */
        else {
            coreblas_error(1, "The case Plasma[Conj]Trans / PlasmaRowwise is not yet implemented");
            return -1;
        }
    }

    return PLASMA_SUCCESS;
}
