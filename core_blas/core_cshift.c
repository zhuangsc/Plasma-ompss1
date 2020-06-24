/**
 *
 * @file core_cshift.c
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
 * @generated c Tue Jan  7 11:44:47 2014
 *
 **/

#include <stdlib.h>
#include "common.h"
#include "quark.h"

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_cshiftw Shift a linear chain of block using a supplied workspace
 *      by following the cycle defined by:  k_(i+1) = (k_i * m) % q;
 *
 *******************************************************************************
 *
 * @param[in] s
 *         Start value in the cycle
 *
 * @param[in] cl
 *         Cycle length
 *         if cl == 0, all the permutations from the cycle are done
 *         else the cycle is split onto several threads and the number
 *         of permutation to do has to be specified to not get overlap
 *
 * @param[in] m
 *         Number of lines of tile A
 *
 * @param[in] n
 *         Number of columns of tile A
 *
 * @param[in] L
 *         Length of each block of data to move
 *
 * @param[in,out] A
 *         Matrix of size m-by-n with each element of size L.
 *         On exit, A = A', where A' contains the permutations
 *
 * @param[in] W
 *         Array of size L. On entry, must contain:
 *         W(:) = A(s*L:s*L+L-1)
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cshiftw = PCORE_cshiftw
#define CORE_cshiftw PCORE_cshiftw
#endif
void CORE_cshiftw(int s, int cl, int m, int n, int L, PLASMA_Complex32_t *A, PLASMA_Complex32_t *W) {
    int64_t k, k1;
    int     i, q, kL, k1L;
    size_t  chunk_size = L*sizeof(PLASMA_Complex32_t);

    q = m * n - 1;
    k = s;

    if( cl != 0 ) {
        for (i=1; i<cl; i++) {
            k1 = (k * m) % (int64_t)q;

            /* A(k*L:k*L+L-1) = A(k1*L:k1*L+L-1) */
            kL  = k *L;
            k1L = k1*L;

            memcpy( A+kL, A+k1L, chunk_size);
            k = k1;
        }
    }
    else {
        while (1) {
            k1 = (k * m) % (int64_t)q;
            if( k1 == s )
                break;

            /* A(k*L:k*L+L-1) = A(k1*L:k1*L+L-1) */
            kL  = k *L;
            k1L = k1*L;
            memcpy( A+kL, A+k1L, chunk_size);
            k = k1;
        }
    }
    kL = k*L;
    memcpy(A+kL, W, chunk_size);
}

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 *  CORE_cshift Shift a cycle of block. Same as core_cshiftw but you
 *     don't need to provide the workspace.  As a matter of fact, the
 *     cycle cannot be split anymore to keep data coherency.
 *
 *******************************************************************************
 *
 * @param[in] s
 *         Start value in the cycle
 *
 * @param[in] m
 *         Number of lines of tile A
 *
 * @param[in] n
 *         Number of columns of tile A
 *
 * @param[in] L
 *         Length of each block of data to move
 *
 * @param[in,out] A
 *         Matrix of size m-by-n with each element of size L.
 *         On exit, A = A', where A' contains the permutations
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cshift = PCORE_cshift
#define CORE_cshift PCORE_cshift
#endif
void CORE_cshift(int s, int m, int n, int L, PLASMA_Complex32_t *A) {
    PLASMA_Complex32_t *W;

    W = (PLASMA_Complex32_t*)malloc(L * sizeof(PLASMA_Complex32_t));
    memcpy(W, &(A[s*L]), L*sizeof(PLASMA_Complex32_t));
    CORE_cshiftw(s, 0, m, n, L, A, W);
    free(W);
}
