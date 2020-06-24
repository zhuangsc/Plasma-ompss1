/**
 *
 * @file core_dgeqp3_pivot.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:44:49 2014
 *
 **/
#include <math.h>
#include <cblas.h>
#include "common.h"

#define A(m,n) BLKADDR( A, double, m, n )

/***************************************************************************//**
 *
 * @ingroup CORE_double
 *
 * CORE_dgeqp3_pivot finds next pivot, pvt, based on maximum column norm.
 * It applies the swap to the matrices A, F, and vectors jpvt, norms1, norms2.
 * If info != 0, it returns immediately, doing no work.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *         On entry, descriptor for m by n matrix A.
 *         On exit, column k of jj-th block column is swapped with column pvt.
 *
 * @param[in,out] F
 *         On entry, n by nb matrix F.
 *         On exit, row k is swapped with row pvt - jj*nb.
 *         Currently, F is stored column-wise, not tile-wise.
 *
 * @param[in] ldf
 *         Leading dimension of F. ldf >= max(1,A.n).
 *
 * @param[in] jj
 *         Index of current block column, 0 <= jj < A.nt.
 *
 * @param[in] k
 *         Index of current column within block column, 0 <= k < A.nb.
 *
 * @param[in,out] jpvt
 *         Permutation vector, dimension n.
 *         On exit, swaps entries jpvt[k+jj*nb] and jpvt[pvt].
 *
 * @param[in,out] norms1
 *         On entry, vector of partial column norms, dimension n.
 *         On exit, sets norms1[pvt] = norms1[k+jj*nb].
 *
 * @param[in,out] norms2
 *         On entry, vector of original column norms, dimension n.
 *         On exit, sets norms2[pvt] = norms2[k+jj*nb].
 *
 * @param[in] info
 *         Error code from dgeqp3_update; zero if no error.
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgeqp3_pivot = PCORE_dgeqp3_pivot
#define CORE_dgeqp3_pivot PCORE_dgeqp3_pivot
#endif
void CORE_dgeqp3_pivot( PLASMA_desc A, double *F, int ldf,
                        int jj, int k, int *jpvt,
                        double *norms1, double *norms2, int *info )
{
    double *Aij, *Aip;
    int pvt, ii, pp, p, mb, lda, tmp, jk;

    /* since pivot depends on all the norm updates, check info here
     * to detect errors from dgeqp3_updates. */
    if ( *info != 0 )
        return;

    /* jk and pvt are indices in global matrix */
    jk = jj*A.nb + k;
    pvt = jk + cblas_idamax( A.n - jk, &norms1[jk], 1 );
    if ( pvt != jk ) {
        tmp       = jpvt[jk];
        jpvt[jk]  = jpvt[pvt];
        jpvt[pvt] = tmp;

        norms1[pvt] = norms1[jk];  /* don't need to save norms[pvt] */
        norms2[pvt] = norms2[jk];

        cblas_dswap( A.nb, &F[k], ldf, &F[pvt-jj*A.nb], ldf );  /* rows k and pvt */

        pp = pvt / A.nb;  /* tile containing pvt  */
        p  = pvt % A.nb;  /* index within pp tile */
        for( ii = 0; ii < A.mt; ++ii ) {
            mb  = min( A.mb, A.m - ii*A.mb );
            lda = BLKLDD( A, ii );
            Aij = A(ii,jj);
            Aip = A(ii,pp);
            cblas_dswap( mb, &Aij[k*lda], 1, &Aip[p*lda], 1 );  /* cols k and pvt */
        }
    }
}
