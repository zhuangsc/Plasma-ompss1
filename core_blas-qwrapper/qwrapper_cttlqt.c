/**
 *
 * @file qwrapper_cttlqt.c
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
void QUARK_CORE_cttlqt(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int ib, int nb,
                       PLASMA_Complex32_t *A1, int lda1,
                       PLASMA_Complex32_t *A2, int lda2,
                       PLASMA_Complex32_t *T, int ldt)
{
    DAG_CORE_TTLQT;
    QUARK_Insert_Task(quark, CORE_cttlqt_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A1,            INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A2,            INOUT|QUARK_REGION_D|QUARK_REGION_L|LOCALITY,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(PLASMA_Complex32_t)*ib*nb,    T,             OUTPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex32_t)*nb,       NULL,          SCRATCH,
        sizeof(PLASMA_Complex32_t)*ib*nb,    NULL,          SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cttlqt_quark = PCORE_cttlqt_quark
#define CORE_cttlqt_quark PCORE_cttlqt_quark
#endif
void CORE_cttlqt_quark(Quark *quark)
{
    int m;
    int n;
    int ib;
    PLASMA_Complex32_t *A1;
    int lda1;
    PLASMA_Complex32_t *A2;
    int lda2;
    PLASMA_Complex32_t *T;
    int ldt;
    PLASMA_Complex32_t *TAU;
    PLASMA_Complex32_t *WORK;

    quark_unpack_args_11(quark, m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
    CORE_cttlqt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
}
