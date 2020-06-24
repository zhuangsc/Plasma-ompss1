/**
 *
 * @file qwrapper_dtsmqr_sytra1.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:44:59 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dtsmqr_sytra1(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_enum side, PLASMA_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              double *A1, int lda1,
                              double *A2, int lda2,
                              const double *V, int ldv,
                              const double *T, int ldt)
{
    int ldwork = side == PlasmaLeft ? ib : nb;

    DAG_CORE_TSMQR;
    QUARK_Insert_Task(quark, CORE_dtsmqr_sytra1_quark, task_flags,
        sizeof(PLASMA_enum),                &side,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(double)*nb*nb,    A1,            INOUT|QUARK_REGION_L|QUARK_REGION_D,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(double)*nb*nb,    A2,            INOUT,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(double)*nb*nb,    V,             INPUT,
        sizeof(int),                        &ldv,   VALUE,
        sizeof(double)*ib*nb,    T,             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(double)*ib*nb,    NULL,          SCRATCH,
        sizeof(int),                        &ldwork, VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dtsmqr_sytra1_quark = PCORE_dtsmqr_sytra1_quark
#define CORE_dtsmqr_sytra1_quark PCORE_dtsmqr_sytra1_quark
#endif
void CORE_dtsmqr_sytra1_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    double *A1;
    int lda1;
    double *A2;
    int lda2;
    double *V;
    int ldv;
    double *T;
    int ldt;
    double *WORK;
    int ldwork;

    quark_unpack_args_18(quark, side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
    CORE_dtsmqr_sytra1(side, trans, m1, n1, m2, n2, k, ib, A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}
