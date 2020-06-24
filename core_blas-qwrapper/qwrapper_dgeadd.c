/**
 *
 * @file qwrapper_dgeadd.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:44:56 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgeadd(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb, double alpha,
                       const double *A, int lda,
                             double *B, int ldb)
{
    DAG_CORE_GEADD;
    QUARK_Insert_Task(quark, CORE_dgeadd_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(double),         &alpha, VALUE,
        sizeof(double)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(double)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgeadd_quark = PCORE_dgeadd_quark
#define CORE_dgeadd_quark PCORE_dgeadd_quark
#endif
void CORE_dgeadd_quark(Quark *quark)
{
    int M;
    int N;
    double alpha;
    double *A;
    int LDA;
    double *B;
    int LDB;

    quark_unpack_args_7(quark, M, N, alpha, A, LDA, B, LDB);
    CORE_dgeadd(M, N, alpha, A, LDA, B, LDB);
    return;
}

