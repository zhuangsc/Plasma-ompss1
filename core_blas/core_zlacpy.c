/**
 *
 * @file core_zlacpy.c
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
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_PLASMA_zlacpycopies all or part of a two-dimensional matrix A to another
 *  matrix B
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = PlasmaUpperLower: All the matrix A
 *            = PlasmaUpper: Upper triangular part
 *            = PlasmaLower: Lower triangular part
 *
 * @param[in] M
 *          The number of rows of the matrices A and B. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrices A and B. N >= 0.
 *
 * @param[in] A
 *          The M-by-N matrix to copy.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] B
 *          The M-by-N copy of the matrix A.
 *          On exit, B = A ONLY in the locations specified by uplo.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlacpy = PCORE_zlacpy
#define CORE_zlacpy PCORE_zlacpy
#endif
void CORE_zlacpy(PLASMA_enum uplo, int M, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 PLASMA_Complex64_t *B, int LDB)
{
    LAPACKE_zlacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}
