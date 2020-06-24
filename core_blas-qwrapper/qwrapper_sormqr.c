/**
 *
 * @file qwrapper_sormqr.c
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
 * @generated s Tue Jan  7 11:44:57 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_sormqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const float *A, int lda,
                       const float *T, int ldt,
                       float *C, int ldc)
{
    DAG_CORE_UNMQR;
    QUARK_Insert_Task(quark, CORE_sormqr_quark, task_flags,
        sizeof(PLASMA_enum),                &side,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(float)*nb*nb,   A,      INPUT | QUARK_REGION_L,
        sizeof(int),                        &lda,   VALUE,
        sizeof(float)*ib*nb,   T,      INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(float)*nb*nb,   C,      INOUT,
        sizeof(int),                        &ldc,   VALUE,
        sizeof(float)*ib*nb,   NULL,   SCRATCH,
        sizeof(int),                        &nb,    VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_sormqr_quark = PCORE_sormqr_quark
#define CORE_sormqr_quark PCORE_sormqr_quark
#endif
void CORE_sormqr_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum trans;
    int m;
    int n;
    int k;
    int ib;
    float *A;
    int lda;
    float *T;
    int ldt;
    float *C;
    int ldc;
    float *WORK;
    int ldwork;

    quark_unpack_args_14(quark, side, trans, m, n, k, ib,
                         A, lda, T, ldt, C, ldc, WORK, ldwork);
    CORE_sormqr(side, trans, m, n, k, ib,
                A, lda, T, ldt, C, ldc, WORK, ldwork);
}
