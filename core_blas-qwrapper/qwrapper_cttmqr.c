/**
 *
 * @file qwrapper_cttmqr.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @generated c Tue Jan  7 11:44:55 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_cttmqr(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum side, PLASMA_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       const PLASMA_Complex32_t *V, int ldv,
                       const PLASMA_Complex32_t *T, int ldt)
{
    int ldwork = side == PlasmaLeft ? ib : nb;

    DAG_CORE_TTMQR;
    QUARK_Insert_Task(quark, CORE_cttmqr_quark, task_flags,
        sizeof(PLASMA_enum),                &side,  VALUE,
        sizeof(PLASMA_enum),                &trans, VALUE,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A1,            INOUT,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A2,            INOUT,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    V,             INPUT|QUARK_REGION_D|QUARK_REGION_U,
        sizeof(int),                        &ldv,   VALUE,
        sizeof(PLASMA_Complex32_t)*ib*nb,    T,             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex32_t)*ib*nb,    NULL,          SCRATCH,
        sizeof(int),                        &ldwork,    VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cttmqr_quark = PCORE_cttmqr_quark
#define CORE_cttmqr_quark PCORE_cttmqr_quark
#endif
void CORE_cttmqr_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    PLASMA_Complex32_t *A1;
    int lda1;
    PLASMA_Complex32_t *A2;
    int lda2;
    PLASMA_Complex32_t *V;
    int ldv;
    PLASMA_Complex32_t *T;
    int ldt;
    PLASMA_Complex32_t *WORK;
    int ldwork;

    quark_unpack_args_18(quark, side, trans, m1, n1, m2, n2, k, ib,
                         A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
    CORE_cttmqr(side, trans, m1, n1, m2, n2, k, ib,
                A1, lda1, A2, lda2, V, ldv, T, ldt, WORK, ldwork);
}
