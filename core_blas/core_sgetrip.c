/**
 *
 * @file core_sgetrip.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @generated s Tue Jan  7 11:44:47 2014
 *
 **/
#include "common.h"

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 * @ingroup CORE_float
 *
 *  CORE_sgetrip transposes a m-by-n matrix in place using an extra
 *      workspace of size m-by-n.
 *      Note : For square tile, workspace is not used.
 *
 *******************************************************************************
 *
 * @param[in] m
 *         Number of lines of tile A
 *
 * @param[in] n
 *         Number of columns of tile A
 *
 * @param[in,out] A
 *         Tile of size m-by-n
 *         On exit, A = trans(A)
 *
 * @param[out] W
 *         Workspace of size n-by-m if n != m, NULL otherwise.
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgetrip = PCORE_sgetrip
#define CORE_sgetrip PCORE_sgetrip
#endif
void CORE_sgetrip(int m, int n, float *A, float *W) {
    float t;
    int i, j;

    if( m != n ) {
        /* rectangular transposition (use workspace) */
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                W[j+i*n] = A[i+j*m];
            }
        }
        memcpy(A, W, m*n*sizeof(float));
    }
    else {
        /* square transposition (swap pairwise) */
        for (i=0; i<m; i++) {
            for (j=i+1; j<n; j++) {
                t        = A[j+i*n];
                A[j+i*n] = A[i+j*m];
                A[i+j*m] = t;
            }
        }
    }
}
