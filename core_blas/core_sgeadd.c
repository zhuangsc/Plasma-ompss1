/**
 *
 * @file core_sgeadd.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
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
 *  CORE_sgeadd adds to matrices together.
 *
 *       B <- alpha * A  + B
 *
 *******************************************************************************
 *
 * @param[in] M
 *          Number of rows of the matrices A and B.
 *
 * @param[in] N
 *          Number of columns of the matrices A and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,M)
 *
 * @param[in,out] B
 *          Matrix of size LDB-by-N.
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgeadd = PCORE_sgeadd
#define CORE_sgeadd PCORE_sgeadd
#endif
int CORE_sgeadd(int M, int N, float alpha,
                const float *A, int LDA,
                      float *B, int LDB)
{
    int j;

    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if ( (LDA < max(1,M)) && (M > 0) ) {
        coreblas_error(5, "Illegal value of LDA");
        return -5;
    }
    if ( (LDB < max(1,M)) && (M > 0) ) {
        coreblas_error(7, "Illegal value of LDB");
        return -7;
    }

    if (M == LDA && M == LDB)
        cblas_saxpy(M*N, (alpha), A, 1, B, 1);
    else {
        for (j = 0; j < N; j++)
            cblas_saxpy(M, (alpha), A + j*LDA, 1, B + j*LDB, 1);
    }

    return PLASMA_SUCCESS;
}
