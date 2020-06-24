/**
 *
 * @file core_dtrssq.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:44:47 2014
 *
 **/
#include <math.h>
#include <lapacke.h>
#include "common.h"

#define REAL

#define UPDATE( __nb, __value )                                         \
    if (__value != 0. ){                                                \
        if ( *scale < __value ) {                                       \
            *sumsq = __nb + (*sumsq) * ( *scale / __value ) * ( *scale / __value ); \
            *scale = __value;                                           \
        } else {                                                        \
            *sumsq = *sumsq + __nb * ( __value / *scale ) *  ( __value / *scale ); \
        }                                                               \
    }

/*****************************************************************************
 *
 * @ingroup CORE_double
 *
 *  CORE_dtrssq returns the values scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( i, j )**2 ) + ( scale**2 )*sumsq,
 *                     i,j
 *
 * with i from 0 to M-1 and j form 0 to N-1. The value of sumsq is
 * assumed to be at least unity and the value of ssq will then satisfy
 *
 *    1.0 .le. ssq .le. ( sumsq + 2*m*n ).
 *
 * scale is assumed to be non-negative and scl returns the value
 *
 *    scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
 *           i
 *
 * scale and sumsq must be supplied in SCALE and SUMSQ respectively.
 * SCALE and SUMSQ are overwritten by scl and ssq respectively.
 *
 * The routine makes only one pass through the tile A.
 * See also LAPACK dlassq.f
 *
 *******************************************************************************
 *
 *  @param[in] M
 *          The number of rows in the tile A.
 *
 *  @param[in] N
 *          The number of columns in the tile A.
 *
 *  @param[in] A
 *          The M-by-N matrix on which to compute the norm.
 *
 *  @param[in] LDA
 *          The leading dimension of the tile A. LDA >= max(1,M).
 *
 *  @param[in,out] scale
 *          On entry, the value  scale  in the equation above.
 *          On exit, scale is overwritten with the value scl.
 *
 *  @param[in,out] sumsq
 *          On entry, the value  sumsq  in the equation above.
 *          On exit, SUMSQ is overwritten with the value ssq.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval -k, the k-th argument had an illegal value
 *
 */
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtrssq = PCORE_dtrssq
#define CORE_dtrssq PCORE_dtrssq
#endif
int CORE_dtrssq(PLASMA_enum uplo, PLASMA_enum diag, int M, int N,
                const double *A, int LDA,
                double *scale, double *sumsq)
{
    int i, j, imax;
    int idiag = (diag == PlasmaUnit) ? 1 : 0;
    double tmp;
    double *ptr;

    if ( diag == PlasmaUnit ){
        tmp = sqrt( min(M, N) );
        UPDATE( 1., tmp );
    }

    if  (uplo == PlasmaUpper ) {
        M = min(M, N);

        for(j=0; j<N; j++) {
            ptr = (double*) ( A + j * LDA );
            imax = min(j+1-idiag, M);

            for(i=0; i<imax; i++, ptr++) {
                tmp = fabs(*ptr);
                UPDATE( 1., tmp );

#ifdef COMPLEX
                ptr++;
                tmp = fabs(*ptr);
                UPDATE( 1., tmp );
#endif
            }
        }
    }
    else {
        N = min(M, N);

        for(j=0; j<N; j++) {
            ptr = (double*) ( A + j * (LDA+1) + idiag );

            for(i=j+idiag; i<M; i++, ptr++) {
                tmp = fabs(*ptr);
                UPDATE( 1., tmp );

#ifdef COMPLEX
                ptr++;
                tmp = fabs(*ptr);
                UPDATE( 1., tmp );
#endif
            }
        }
    }
    return PLASMA_SUCCESS;
}
