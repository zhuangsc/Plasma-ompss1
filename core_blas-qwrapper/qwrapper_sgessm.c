/**
 *
 * @file qwrapper_sgessm.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:44:56 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sgessm(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int k, int ib, int nb,
                       const int *IPIV,
                       const float *L, int ldl,
                       float *A, int lda)
{
    DAG_CORE_GESSM;
    QUARK_Insert_Task(quark, CORE_sgessm_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(int)*nb,                      IPIV,          INPUT,
        sizeof(float)*nb*nb,    L,             INPUT | QUARK_REGION_L,
        sizeof(int),                        &ldl,   VALUE,
        sizeof(float)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sgessm_quark = PCORE_sgessm_quark
#define CORE_sgessm_quark PCORE_sgessm_quark
#endif
void CORE_sgessm_quark(Quark *quark)
{
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    float *L;
    int ldl;
    float *A;
    int lda;

    quark_unpack_args_9(quark, m, n, k, ib, IPIV, L, ldl, A, lda);
    CORE_sgessm(m, n, k, ib, IPIV, L, ldl, A, lda);
}
