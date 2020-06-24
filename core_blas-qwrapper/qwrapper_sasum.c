/**
 *
 * @file qwrapper_sasum.c
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

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                       const float *A, int lda, int szeA,
                       float *work, int szeW)
{
    QUARK_Insert_Task(
        quark, CORE_sasum_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(float)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*szeW,                 work,              INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sasum_quark = PCORE_sasum_quark
#define CORE_sasum_quark PCORE_sasum_quark
#endif
void CORE_sasum_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    int M;
    int N;
    float *A;
    int lda;
    float *work;

    quark_unpack_args_7(quark, storev, uplo, M, N, A, lda, work);
    CORE_sasum(storev, uplo, M, N, A, lda, work);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                          const float *A, int lda, int szeA,
                          float *work, int szeW, float *fake, int szeF)
{
    DAG_CORE_ASUM;
    if ( work == fake ) {
        QUARK_Insert_Task(
            quark, CORE_sasum_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(float)*szeA,     A,                 INPUT,
            sizeof(int),                        &lda,       VALUE,
            sizeof(float)*szeW,                 work,              INOUT | GATHERV,
            sizeof(float)*1,                    NULL,              SCRATCH,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_sasum_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(float)*szeA,     A,                 INPUT,
            sizeof(int),                        &lda,       VALUE,
            sizeof(float)*szeW,                 work,              INOUT,
            sizeof(float)*szeF,                 fake,              OUTPUT | GATHERV,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sasum_f1_quark = PCORE_sasum_f1_quark
#define CORE_sasum_f1_quark PCORE_sasum_f1_quark
#endif
void CORE_sasum_f1_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    int M;
    int N;
    float *A;
    int lda;
    float *work;
    float *fake;

    quark_unpack_args_8(quark, storev, uplo, M, N, A, lda, work, fake);
    CORE_sasum(storev, uplo, M, N, A, lda, work);
}
