/**
 *
 * @file qwrapper_cgeadd.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Tue Jan  7 11:44:56 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cgeadd(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb, PLASMA_Complex32_t alpha,
                       const PLASMA_Complex32_t *A, int lda,
                             PLASMA_Complex32_t *B, int ldb)
{
    DAG_CORE_GEADD;
    QUARK_Insert_Task(quark, CORE_cgeadd_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cgeadd_quark = PCORE_cgeadd_quark
#define CORE_cgeadd_quark PCORE_cgeadd_quark
#endif
void CORE_cgeadd_quark(Quark *quark)
{
    int M;
    int N;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int LDA;
    PLASMA_Complex32_t *B;
    int LDB;

    quark_unpack_args_7(quark, M, N, alpha, A, LDA, B, LDB);
    CORE_cgeadd(M, N, alpha, A, LDA, B, LDB);
    return;
}

