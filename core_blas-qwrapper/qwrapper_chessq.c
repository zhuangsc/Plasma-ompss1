/**
 *
 * @file qwrapper_chessq.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Tue Jan  7 11:44:57 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_chessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum uplo, int n, const PLASMA_Complex32_t *A, int lda,
                           float *scale, float *sumsq,
                           float *fake, int szeF, int paramF )
{
    DAG_CORE_LASSQ;

    if ( (fake == scale) && (paramF & GATHERV) ) {
        QUARK_Insert_Task(quark, CORE_chessq_quark, task_flags,
            sizeof(PLASMA_enum),              &uplo, VALUE,
            sizeof(int),                      &n,    VALUE,
            sizeof(PLASMA_Complex32_t)*lda*n, A,         INPUT,
            sizeof(int),                      &lda,  VALUE,
            sizeof(float)*1,                 scale,     INOUT | paramF,
            sizeof(float)*1,                 sumsq,     INOUT,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_chessq_f1_quark, task_flags,
            sizeof(PLASMA_enum),              &uplo, VALUE,
            sizeof(int),                      &n,    VALUE,
            sizeof(PLASMA_Complex32_t)*lda*n, A,         INPUT,
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
#pragma weak CORE_chessq_quark = PCORE_chessq_quark
#define CORE_chessq_quark PCORE_chessq_quark
#endif
void CORE_chessq_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int n;
    PLASMA_Complex32_t *A;
    int lda;
    float *scale;
    float *sumsq;

    quark_unpack_args_6( quark, uplo, n, A, lda, scale, sumsq );
    CORE_chessq( uplo, n, A, lda, scale, sumsq );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_chessq_f1_quark = PCORE_chessq_f1_quark
#define CORE_chessq_f1_quark PCORE_chessq_f1_quark
#endif
void CORE_chessq_f1_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int n;
    PLASMA_Complex32_t *A;
    int lda;
    float *scale;
    float *sumsq;
    float *fake;

    quark_unpack_args_7( quark, uplo, n, A, lda, scale, sumsq, fake );
    CORE_chessq( uplo, n, A, lda, scale, sumsq );
}
