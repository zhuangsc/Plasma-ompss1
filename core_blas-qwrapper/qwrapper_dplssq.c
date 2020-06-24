/**
 *
 * @file qwrapper_dplssq.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:44:57 2014
 *
 **/
#include <math.h>
#include "common.h"

/*****************************************************************************
 *
 * @ingroup CORE_double
 *
 *  QUARK_CORE_dplssq returns: scl * sqrt(ssq)
 *
 * with scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( 2*i )**2 * A( 2*i+1 ) )
 *                      i
 *
 * The values of A(2*i+1) are assumed to be at least unity.
 * The values of A(2*i) are assumed to be non-negative and scl is
 *
 *    scl = max( A( 2*i ) ),
 *           i
 *
 * The routine makes only one pass through the matrix A.
 *
 *******************************************************************************
 *
 *  @param[in] M
 *          The number of couple (scale, sumsq) in the matrix A.
 *
 *  @param[in] A
 *          The 2-by-M matrix.
 *
 *  @param[out] result
 *          On exit, result contains scl * sqrt( ssq )
 *
 */
void QUARK_CORE_dplssq( Quark *quark, Quark_Task_Flags *task_flags,
                        int m, const double *A, double *result )
{
    DAG_CORE_LASSQ;
    QUARK_Insert_Task(quark, CORE_dplssq_quark, task_flags,
        sizeof(int),       &m,      VALUE,
        sizeof(double)*m,   A,          INOUT,
        sizeof(double)*1,   result,     OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dplssq_quark = PCORE_dplssq_quark
#define CORE_dplssq_quark PCORE_dplssq_quark
#endif
void CORE_dplssq_quark(Quark *quark)
{
    int m, i;
    double *A;
    double *result;

    quark_unpack_args_3( quark, m, A, result );

    for( i=1; i<m; i++ ) {
        if( A[0] < A[i] ) {
            A[1] = A[2*i+1] + (A[1]     * (( A[0] / A[2*i] ) * ( A[0] / A[2*i] )));
            A[0] = A[2*i];
        } else {
            A[1] = A[1]     + (A[2*i+1] * (( A[2*i] / A[0] ) * ( A[2*i] / A[0] )));
        }
    }

    *result = A[0] * sqrt( A[1] );
}
