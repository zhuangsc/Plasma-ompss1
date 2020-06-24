/**
 *
 * @file core_strmm.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:44:46 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_strmm - Computes B = alpha*op( A )*B or B = alpha*B*op( A ).
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = PlasmaLeft:  B = alpha*op( A )*B.
 *          = PlasmaRight: B = alpha*B*op( A ).
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or ugate transposed:
 *          = PlasmaNoTrans:   A is transposed;
 *          = PlasmaTrans:     A is not transposed;
 *          = PlasmaTrans: A is ugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] M
 *          The number of rows of the matrix B. M >= 0.
 *
 * @param[in] N
 *          The number of columns pf the matrix B. N >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is an array of dimansion LDA-by-k, where k = M if side =
 *          PlasmaLeft and k =N when side = PlasmaRight.
 *          The triangular matrix A. If uplo = PlasmaUpper, the leading k-by-k upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = PlasmaLower, the leading k-by-k
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = PlasmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,k). Where k = M if side =
 *          PlasmaLeft and k =N when side = PlasmaRight.
 *
 * @param[in,out] B
 *          On entry, the LDB-by-N matrix B.
 *          On exit, if return value = 0, the M-by-N matrix is overwritten by
 *          the transformed matrix.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_strmm = PCORE_strmm
#define CORE_strmm PCORE_strmm
#endif
void CORE_strmm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                float alpha,
                const float *A, int LDA,
                float *B, int LDB)
{
    cblas_strmm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        (alpha), A, LDA,
        B, LDB);
}
