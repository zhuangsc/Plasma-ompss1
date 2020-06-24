/**
 *
 * @file qwrapper_strssq.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:44:57 2014
 *
 **/
#include "common.h"

void
CORE_strssq_quark(Quark *quark);
void
CORE_strssq_f1_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_strssq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, PLASMA_enum diag,
                           int m, int n, const float *A, int lda,
                           float *scale, float *sumsq,
                           float *fake, int szeF, int paramF )
{
    DAG_CORE_LASSQ;
    if ( (fake == scale) && (paramF & GATHERV) ) {
        QUARK_Insert_Task(quark, CORE_strssq_quark, task_flags,
            sizeof(PLASMA_enum),              &uplo, VALUE,
            sizeof(PLASMA_enum),              &diag, VALUE,
            sizeof(int),                      &m,    VALUE,
            sizeof(int),                      &n,    VALUE,
            sizeof(float)*lda*n, A,         INPUT,
            sizeof(int),                      &lda,  VALUE,
            sizeof(float)*1,                 scale,     INOUT | GATHERV,
            sizeof(float)*1,                 sumsq,     INOUT,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_strssq_f1_quark, task_flags,
            sizeof(PLASMA_enum),              &uplo, VALUE,
            sizeof(PLASMA_enum),              &diag, VALUE,
            sizeof(int),                      &m,    VALUE,
            sizeof(int),                      &n,    VALUE,
            sizeof(float)*lda*n, A,         INPUT,
            sizeof(int),                      &lda,  VALUE,
            sizeof(float)*1,                 scale,     INOUT,
            sizeof(float)*1,                 sumsq,     INOUT,
            sizeof(float)*szeF,              fake,      paramF,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_strssq_quark = PCORE_strssq_quark
#define CORE_strssq_quark PCORE_strssq_quark
#endif
void CORE_strssq_quark(Quark *quark)
{
    PLASMA_enum uplo, diag;
    int m;
    int n;
    float *A;
    int lda;
    float *scale;
    float *sumsq;

    quark_unpack_args_8( quark, uplo, diag, m, n, A, lda, scale, sumsq );
    CORE_strssq( uplo, diag, m, n, A, lda, scale, sumsq );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_strssq_f1_quark = PCORE_strssq_f1_quark
#define CORE_strssq_f1_quark PCORE_strssq_f1_quark
#endif
void CORE_strssq_f1_quark(Quark *quark)
{
    PLASMA_enum uplo, diag;
    int m;
    int n;
    float *A;
    int lda;
    float *scale;
    float *sumsq;
    float *fake;

    quark_unpack_args_9( quark, uplo, diag, m, n, A, lda, scale, sumsq, fake );
    CORE_strssq( uplo, diag, m, n, A, lda, scale, sumsq );
}
