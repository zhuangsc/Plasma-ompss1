/**
 *
 * @file qwrapper_clacpy.c
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
void QUARK_CORE_clacpy(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int m, int n, int nb,
                       const PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex32_t *B, int ldb)
{
    DAG_CORE_LACPY;
    QUARK_Insert_Task(quark, CORE_clacpy_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,             OUTPUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clacpy_quark = PCORE_clacpy_quark
#define CORE_clacpy_quark PCORE_clacpy_quark
#endif
void CORE_clacpy_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    const PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t *B;
    int LDB;

    quark_unpack_args_7(quark, uplo, M, N, A, LDA, B, LDB);
    LAPACKE_clacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clacpy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum uplo, int m, int n, int nb,
                          const PLASMA_Complex32_t *A, int lda,
                          PLASMA_Complex32_t *B, int ldb,
                          PLASMA_Complex32_t *fake1, int szefake1, int flag1)
{
    DAG_CORE_LACPY;
    if ( fake1 == B ) {
        QUARK_Insert_Task(quark, CORE_clacpy_quark, task_flags,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(int),                        &m,     VALUE,
            sizeof(int),                        &n,     VALUE,
            sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT,
            sizeof(int),                        &lda,   VALUE,
            sizeof(PLASMA_Complex32_t)*nb*nb,    B,             OUTPUT | flag1,
            sizeof(int),                        &ldb,   VALUE,
            0);
    }
    else {
        QUARK_Insert_Task(quark, CORE_clacpy_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(int),                        &m,     VALUE,
            sizeof(int),                        &n,     VALUE,
            sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT,
            sizeof(int),                        &lda,   VALUE,
            sizeof(PLASMA_Complex32_t)*nb*nb,    B,             OUTPUT,
            sizeof(int),                        &ldb,   VALUE,
            sizeof(PLASMA_Complex32_t)*szefake1, fake1,         flag1,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clacpy_f1_quark = PCORE_clacpy_f1_quark
#define CORE_clacpy_f1_quark PCORE_clacpy_f1_quark
#endif
void CORE_clacpy_f1_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    const PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t *B;
    int LDB;
    void *fake1;

    quark_unpack_args_8(quark, uplo, M, N, A, LDA, B, LDB, fake1);
    LAPACKE_clacpy_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, A, LDA, B, LDB);
}

