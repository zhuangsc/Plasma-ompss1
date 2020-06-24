/**
 *
 * @file core_spltmg_circul.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:44:48 2014
 *
 **/
#include <lapacke.h>
#include <math.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_spltmg_circul is a kernel used in circulant matrix generation
 *
 *  See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999880
 *
 *  Circulant matrix
 *
 *  A circulant matrix has the property that each row is obtained from the
 *  previous one by cyclically permuting the entries one step forward. It is a
 *  special Toeplitz matrix in which the diagonals "wrap around."
 *
 *  The eigensystem of C (n-by-n) is known explicitly: If t is an nth root of
 *  unity, then the inner product of v and w = [1 t t2 ... t(n – 1)] is an
 *  eigenvalue of C and w(n:-1:1) is an eigenvector, where v is the first column of
 *  C.
 *
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A to initialize. M >= 2.
 *
 * @param[in] N
 *         The number of columns of the tile A to initialize. N >= 0.
 *
 * @param[out] A
 *         On entry, the M-by-N tile to be initialized.
 *
 * @param[in] LDA
 *         The leading dimension of the tile A. LDA >= max(1,M).
 *
 * @param[in] gM
 *         The global number of rows of the full matrix, A is belonging to. gM >= (m0+gM).
 *
 * @param[in] m0
 *         The index of the first row of tile A in the full matrix. m0 >= 0.
 *
 * @param[in] n0
 *         The index of the first column of tile A in the full matrix. n0 >= 0.
 *
 * @param[in] V
 *          Workspace of size gM, that contains the first clumn of the full matrix
 *
 *******************************************************************************
 *
 * @return
 *         \retval PLASMA_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_spltmg_circul = PCORE_spltmg_circul
#define CORE_spltmg_circul PCORE_spltmg_circul
#endif
int CORE_spltmg_circul( int M, int N, float *A, int LDA,
                        int gM, int m0, int n0,
                        const float *V )
{
    int i, j, ii, jj;

    /* Check input arguments */
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if ((LDA < max(1,M)) && (M > 0)) {
        coreblas_error(4, "Illegal value of LDA");
        return -4;
    }
    if (m0 < 0) {
        coreblas_error(6, "Illegal value of m0");
        return -6;
    }
    if (n0 < 0) {
        coreblas_error(7, "Illegal value of n0");
        return -7;
    }
    if (gM < m0+M) {
        coreblas_error(5, "Illegal value of gM");
        return -5;
    }

    for (j=0, jj=n0; j<N; j++, jj++) {
        for (i=0, ii=m0; i<M; i++, ii++) {
            A[LDA*j + i] = V[ (jj-ii+gM)%gM ];
        }
    }

    return PLASMA_SUCCESS;
}
