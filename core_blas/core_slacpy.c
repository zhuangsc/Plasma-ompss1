/**
 *
 * @file core_slacpy.c
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
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_float
 *
 *  CORE_PLASMA_slacpycopies all or part of a two-dimensional matrix A to another
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
#pragma weak CORE_slacpy = PCORE_slacpy
#define CORE_slacpy PCORE_slacpy
#endif
void CORE_slacpy(PLASMA_enum uplo, int M, int N,
                 const float *A, int LDA,
                 float *B, int LDB)
{
    LAPACKE_slacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}
