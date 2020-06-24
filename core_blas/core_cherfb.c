/**
 *
 * @file core_cherfb.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Hatem Ltaief
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
 *  CORE_cherfb overwrites the symmetric complex N-by-N tile C with
 *
 *    Q**T*C*Q
 *
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k)
 *
 *  as returned by CORE_cgeqrt. Only PlasmaLower supported!
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *         @arg PlasmaLower : the upper part of the symmetric matrix C
 *                            is not referenced.
 *         @arg PlasmaUpper : the lower part of the symmetric matrix C
 *                            is not referenced (not supported).
 *
 * @param[in] n
 *         The number of rows/columns of the tile C.  N >= 0.
 *
 * @param[in] k
 *         The number of elementary reflectors whose product defines
 *         the matrix Q. K >= 0.
 *
 * @param[in] ib
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] nb
 *         The blocking size.  NB >= 0.
 *
 * @param[in] A
 *         The i-th column must contain the vector which defines the
 *         elementary reflector H(i), for i = 1,2,...,k, as returned by
 *         CORE_cgeqrt in the first k columns of its array argument A.
 *
 * @param[in] lda
 *         The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] T
 *         The IB-by-K triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] ldt
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[in,out] C
 *         On entry, the symmetric N-by-N tile C.
 *         On exit, C is overwritten by Q**T*C*Q.
 *
 * @param[in] ldc
 *         The leading dimension of the array C. LDC >= max(1,M).
 *
 * @param[in,out] WORK
 *         On exit, if INFO = 0, WORK(1) returns the optimal LDWORK.
 *
 * @param[in] ldwork
 *         The dimension of the array WORK. LDWORK >= max(1,N);
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cherfb = PCORE_cherfb
#define CORE_cherfb PCORE_cherfb
#define CORE_cunmlq PCORE_cunmlq
#define CORE_cunmqr PCORE_cunmqr
int  CORE_cunmlq(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int IB, int K,
                 const PLASMA_Complex32_t *V, int LDV,
                 const PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *C, int LDC,
                 PLASMA_Complex32_t *WORK, int LDWORK);
int  CORE_cunmqr(PLASMA_enum side, PLASMA_enum trans,
                 int M, int N, int K, int IB,
                 const PLASMA_Complex32_t *V, int LDV,
                 const PLASMA_Complex32_t *T, int LDT,
                 PLASMA_Complex32_t *C, int LDC,
                 PLASMA_Complex32_t *WORK, int LDWORK);
#endif
int CORE_cherfb( PLASMA_enum uplo, int n,
                 int k, int ib, int nb,
                 const PLASMA_Complex32_t *A, int lda,
                 const PLASMA_Complex32_t *T, int ldt,
                 PLASMA_Complex32_t *C, int ldc,
                 PLASMA_Complex32_t *WORK, int ldwork )
{
    PLASMA_Complex32_t tmp;
    int i, j;

    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        coreblas_error(1, "Illegal value of uplo");
        return -1;
    }
    if (n < 0) {
        coreblas_error(2, "Illegal value of n");
        return -2;
    }
    if (k < 0) {
        coreblas_error(3, "Illegal value of k");
        return -3;
    }
    if (ib < 0) {
        coreblas_error(4, "Illegal value of ib");
        return -4;
    }
    if (nb < 0) {
        coreblas_error(5, "Illegal value of nb");
        return -5;
    }
    if ( (lda < max(1,n)) && (n > 0) ) {
        coreblas_error(7, "Illegal value of lda");
        return -7;
    }
    if ( (ldt < max(1,ib)) && (ib > 0) ) {
        coreblas_error(9, "Illegal value of ldt");
        return -9;
    }
    if ( (ldc < max(1,n)) && (n > 0) ) {
        coreblas_error(11, "Illegal value of ldc");
        return -11;
    }

    if (uplo == PlasmaLower) {
        /* Rebuild the symmetric block: WORK <- C */
        for (j = 0; j < n; j++) {
            *(WORK + j + j * ldwork) =  *(C + ldc*j + j);
            for (i = j+1; i < n; i++){
                tmp = *(C + i + j*ldc);
                *(WORK + i + j * ldwork) = tmp;
                *(WORK + j + i * ldwork) = conjf( tmp );
            }
        }

        /* Left */
        CORE_cunmqr(PlasmaLeft, PlasmaConjTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        /* Right */
        CORE_cunmqr(PlasmaRight, PlasmaNoTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);

        /*
         * Copy back the final result to the lower part of C
         */
        LAPACKE_clacpy_work( LAPACK_COL_MAJOR, lapack_const(PlasmaLower), n, n, WORK, ldwork, C, ldc );
    }
    else {
        /* Rebuild the symmetric block: WORK <- C */
        for (j = 0; j < n; j++) {
            for (i = 0; i < j; i++){
                tmp = *(C + i + j*ldc);
                *(WORK + i + j * ldwork) = tmp;
                *(WORK + j + i * ldwork) = conjf( tmp );
            }
            *(WORK + j + j * ldwork) =  *(C + ldc*j + j);
        }

        /* Right */
        CORE_cunmlq(PlasmaRight, PlasmaConjTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);
        /* Left */
        CORE_cunmlq(PlasmaLeft, PlasmaNoTrans, n, n, k, ib,
                    A, lda, T, ldt, WORK, ldwork, WORK+nb*ldwork, ldwork);

        /*
         * Copy back the final result to the upper part of C
         */
        LAPACKE_clacpy_work( LAPACK_COL_MAJOR, lapack_const(PlasmaUpper), n, n, WORK, ldwork, C, ldc );
    }
    return 0;
}
