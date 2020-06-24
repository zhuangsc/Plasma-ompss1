/**
 *
 * @file qwrapper_stsmlq_sytra1.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:44:59 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_stsmlq_sytra1(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_enum side, PLASMA_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              float *A1, int lda1,
                              float *A2, int lda2,
                              const float *V, int ldv,
                              const float *T, int ldt)
{
    int ldwork = side == PlasmaLeft ? ib : nb;

    DAG_CORE_TSMLQ;
    QUARK_Insert_Task(quark, CORE_stsmlq_sytra1_quark, task_flags,
        sizeof(PLASMA_enum),                &side,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(float)*nb*nb,    A1,            INOUT|QUARK_REGION_U|QUARK_REGION_D,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(float)*nb*nb,    A2,            INOUT,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(float)*nb*nb,    V,             INPUT,
        sizeof(int),                        &ldv,   VALUE,
        sizeof(float)*ib*nb,    T,             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(float)*ib*nb,    NULL,          SCRATCH,
        sizeof(int),                        &ldwork, VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_stsmlq_sytra1_quark = PCORE_stsmlq_sytra1_quark
#define CORE_stsmlq_sytra1_quark PCORE_stsmlq_sytra1_quark
#endif
void CORE_stsmlq_sytra1_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    float *A1;
    int lda1;
    float *A2;
    int lda2;
    float *V;
    int ldv;
    float *T;
    int ldt;
    float *WORK;
    int ldwork;

    quark_unpack_args_18(quark, side, trans, m1, n1, m2, n2, k, ib,
                         A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
    CORE_stsmlq_sytra1(side, trans, m1, n1, m2, n2, k, ib,
                       A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}

