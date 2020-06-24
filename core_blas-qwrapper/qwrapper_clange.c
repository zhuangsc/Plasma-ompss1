/**
 *
 * @file qwrapper_clange.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Tue Jan  7 11:44:56 2014
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clange(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int M, int N,
                       const PLASMA_Complex32_t *A, int LDA, int szeA,
                       int szeW, float *result)
{
    szeW = max(1, szeW);
    DAG_CORE_LANGE;
    QUARK_Insert_Task(quark, CORE_clange_quark, task_flags,
        sizeof(PLASMA_enum),                &norm,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(PLASMA_Complex32_t)*szeA,     A,             INPUT,
        sizeof(int),                        &LDA,   VALUE,
        sizeof(float)*szeW,                 NULL,          SCRATCH,
        sizeof(float),                      result,        OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clange_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, int M, int N,
                          const PLASMA_Complex32_t *A, int LDA, int szeA,
                          int szeW, float *result,
                          float *fake, int szeF)
{
    szeW = max(1, szeW);
    DAG_CORE_LANGE;

    if ( result == fake ) {
        QUARK_Insert_Task(quark, CORE_clange_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(int),                        &M,     VALUE,
            sizeof(int),                        &N,     VALUE,
            sizeof(PLASMA_Complex32_t)*szeA,     A,             INPUT,
            sizeof(int),                        &LDA,   VALUE,
            sizeof(float)*szeW,                 NULL,          SCRATCH,
            sizeof(float),                      result,        OUTPUT | GATHERV,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_clange_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(int),                        &M,     VALUE,
            sizeof(int),                        &N,     VALUE,
            sizeof(PLASMA_Complex32_t)*szeA,     A,             INPUT,
            sizeof(int),                        &LDA,   VALUE,
            sizeof(float)*szeW,                 NULL,          SCRATCH,
            sizeof(float),                      result,        OUTPUT,
            sizeof(float)*szeF,                 fake,          OUTPUT | GATHERV,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clange_quark = PCORE_clange_quark
#define CORE_clange_quark PCORE_clange_quark
#endif
void CORE_clange_quark(Quark *quark)
{
    float *normA;
    int norm;
    int M;
    int N;
    PLASMA_Complex32_t *A;
    int LDA;
    float *work;

    quark_unpack_args_7(quark, norm, M, N, A, LDA, work, normA);
    *normA = LAPACKE_clange_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        M, N, A, LDA, work);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clange_f1_quark = PCORE_clange_f1_quark
#define CORE_clange_f1_quark PCORE_clange_f1_quark
#endif
void CORE_clange_f1_quark(Quark *quark)
{
    float *normA;
    int norm;
    int M;
    int N;
    PLASMA_Complex32_t *A;
    int LDA;
    float *work;
    float *fake;

    quark_unpack_args_8(quark, norm, M, N, A, LDA, work, normA, fake);
    *normA = LAPACKE_clange_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        M, N, A, LDA, work);
}

