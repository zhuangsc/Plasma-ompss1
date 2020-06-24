/**
 *
 * @file qwrapper_sgeadd.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:44:56 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgeadd(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb, float alpha,
                       const float *A, int lda,
                             float *B, int ldb)
{
    DAG_CORE_GEADD;
    QUARK_Insert_Task(quark, CORE_sgeadd_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(float),         &alpha, VALUE,
        sizeof(float)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(float)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgeadd_quark = PCORE_sgeadd_quark
#define CORE_sgeadd_quark PCORE_sgeadd_quark
#endif
void CORE_sgeadd_quark(Quark *quark)
{
    int M;
    int N;
    float alpha;
    float *A;
    int LDA;
    float *B;
    int LDB;

    quark_unpack_args_7(quark, M, N, alpha, A, LDA, B, LDB);
    CORE_sgeadd(M, N, alpha, A, LDA, B, LDB);
    return;
}

