/**
 *
 * @file core_csyssq.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Tue Jan  7 11:44:47 2014
 *
 **/
#include <math.h>
#include <lapacke.h>
#include "common.h"

#define COMPLEX

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
 * @ingroup CORE_PLASMA_Complex32_t
 *
 *  CORE_csyssq returns the values scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( i, j )**2 ) + ( scale**2 )*sumsq,
 *                     i,j
 *
 * with i from 0 to N-1 and j form 0 to N-1. The value of sumsq is
 * assumed to be at least unity and the value of ssq will then satisfy
 *
 *    1.0 .le. ssq .le. ( sumsq + 2*n*n ).
 *
 * scale is assumed to be non-negative and scl returns the value
 *
 *    scl = max( scale, abs( real( A( i, j ) ) ), abs( aimagf( A( i, j ) ) ) ),
 *          i,j
 *
 * scale and sumsq must be supplied in SCALE and SUMSQ respectively.
 * SCALE and SUMSQ are overwritten by scl and ssq respectively.
 *
 * The routine makes only one pass through the tile triangular part of the
 * symmetric tile A defined by uplo.
 * See also LAPACK classq.f
 *
 *******************************************************************************
 *
 *  @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the symmetric matrix A is to be referenced as follows:
 *          = PlasmaLower:     Only the lower triangular part of the
 *                             symmetric matrix A is to be referenced.
 *          = PlasmaUpper:     Only the upper triangular part of the
 *                             symmetric matrix A is to be referenced.
 *
 *  @param[in] N
 *          The number of columns and rows in the tile A.
 *
 *  @param[in] A
 *          The N-by-N matrix on which to compute the norm.
 *
 *  @param[in] LDA
 *          The leading dimension of the tile A. LDA >= max(1,N).
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
#pragma weak CORE_csyssq = PCORE_csyssq
#define CORE_csyssq PCORE_csyssq
#endif
int CORE_csyssq(PLASMA_enum uplo, int N,
                const PLASMA_Complex32_t *A, int LDA,
                float *scale, float *sumsq)
{
    int i, j;
    float tmp;
    float *ptr;

    if ( uplo == PlasmaUpper ) {
        for(j=0; j<N; j++) {
            ptr = (float*) ( A + j * LDA );

            for(i=0; i<j; i++, ptr++) {

                tmp = fabs(*ptr);
                UPDATE( 2., tmp );

#ifdef COMPLEX
                ptr++;
                tmp = fabs(*ptr);
                UPDATE( 2., tmp );
#endif
            }

            /* Diagonal */
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );

#ifdef COMPLEX
            ptr++;
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );
#endif
        }
    } else {

        for(j=0; j<N; j++) {
            ptr = (float*) ( A + j * LDA + j);

            /* Diagonal */
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );
            ptr++;

#ifdef COMPLEX
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );
            ptr++;
#endif

            for(i=j+1; i<N; i++, ptr++) {

                tmp = fabs(*ptr);
                UPDATE( 2., tmp );

#ifdef COMPLEX
                ptr++;
                tmp = fabs(*ptr);
                UPDATE( 2., tmp );
#endif
            }
        }
    }
    return PLASMA_SUCCESS;
}
