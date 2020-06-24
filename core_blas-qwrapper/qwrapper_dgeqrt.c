/**
 *
 * @file qwrapper_dgeqrt.c
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
 * @generated d Tue Jan  7 11:44:55 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgeqrt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       double *A, int lda,
                       double *T, int ldt)
{
    DAG_CORE_GEQRT;
    QUARK_Insert_Task(quark, CORE_dgeqrt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(double)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(double)*ib*nb,    T,             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(double)*nb,       NULL,          SCRATCH,
        sizeof(double)*ib*nb,    NULL,          SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgeqrt_quark = PCORE_dgeqrt_quark
#define CORE_dgeqrt_quark PCORE_dgeqrt_quark
#endif
void CORE_dgeqrt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    double *A;
    int lda;
    double *T;
    int ldt;
    double *TAU;
    double *WORK;

    quark_unpack_args_9(quark, m, n, ib, A, lda, T, ldt, TAU, WORK);
    CORE_dgeqrt(m, n, ib, A, lda, T, ldt, TAU, WORK);
}
