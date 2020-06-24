/**
 *
 * @file qwrapper_zgeadd.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgeadd(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb, PLASMA_Complex64_t alpha,
                       const PLASMA_Complex64_t *A, int lda,
                             PLASMA_Complex64_t *B, int ldb)
{
    DAG_CORE_GEADD;
    QUARK_Insert_Task(quark, CORE_zgeadd_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeadd_quark = PCORE_zgeadd_quark
#define CORE_zgeadd_quark PCORE_zgeadd_quark
#endif
void CORE_zgeadd_quark(Quark *quark)
{
    int M;
    int N;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;

    quark_unpack_args_7(quark, M, N, alpha, A, LDA, B, LDB);
    CORE_zgeadd(M, N, alpha, A, LDA, B, LDB);
    return;
}

