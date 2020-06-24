/**
 *
 * @file core_sgeqp3_tntpiv.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:44:50 2014
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************
 *
 * @ingroup CORE_float
 *
 *  CORE_sgeqp3_tntpiv computes a QR factorization with column pivoting of a
 *  matrix A:  A*P = Q*R  using Level 3 BLAS.
 *
 *  The matrix Q is represented as a product of elementary reflectors
 *
 *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
 *
 *  Each H(i) has the form
 *
 *     H(i) = I - tau * v * v**T
 *
 *  where tau is a complex scalar, and v is a real/complex vector
 *  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
 *  A(i+1:m,i), and tau in TAU(i).
 *
 *******************************************************************************
 *
 *  Arguments:
 *  ==========
 *
 * @param[in] m
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[in,out] A
 *          A is COMPLEX*16 array, dimension (LDA,N)
 *          On entry, the M-by-N matrix A.
 *          On exit, the upper triangle of the array contains the
 *          min(M,N)-by-N upper trapezoidal matrix R; the elements below
 *          the diagonal, together with the array TAU, represent the
 *          unitary matrix Q as a product of min(M,N) elementary
 *          reflectors.
 *
 * @param[in] lda
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] IPIV
 *          IPIV is INTEGER array, dimension min(M,N)
 *          The pivot indices; for 1 <= j <= min(M,N), column j of the
 *          tile was interchanged with column IPIV(j).
 *
 * @param[out] TAU
 *          TAU is COMPLEX*16 array, dimension (min(M,N))
 *          The scalar factors of the elementary reflectors.
 *
 * @param[in,out] iwork
 *          iwork is INTEGER array, dimension (N)
 *          On entry, if iwork(J).ne.0, the J-th column of A is permuted
 *          to the front of A*P (a leading column); if iwork(J)=0,
 *          the J-th column of A is a free column.
 *          On exit, if iwork(J)=K, then the J-th column of A*P was the
 *          the K-th column of A.
 *
 * @param[out] INFO
 *          = 0: successful exit.
 *          < 0: if INFO = -i, the i-th argument had an illegal value.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgeqp3_tntpiv = PCORE_sgeqp3_tntpiv
#define CORE_sgeqp3_tntpiv PCORE_sgeqp3_tntpiv
#endif
int CORE_sgeqp3_tntpiv(int m, int n,
                       float *A, int lda,
                       int *IPIV, float *tau,
                       int *iwork)
{
    int i, tmp, info;
    memset(iwork, 0, n*sizeof(int));
    info = LAPACKE_sgeqp3(LAPACK_COL_MAJOR, m, n, A, lda, iwork, tau );

    /* Convert IPIV from permutation array, to pivot array
     * WARNING: this is because this kernel is only used in
     * tournament pivoting with rank revealing QR */
    if (info == 0) {
        for(i=0; i<min(m,n); i++) {
            assert(iwork[i] != 0 );

            tmp = iwork[i]-1;
            while( tmp < i ) {
                tmp = IPIV[ tmp ] - 1;
            }
            IPIV[i] = tmp+1;
        }
    }
    return info;
}
