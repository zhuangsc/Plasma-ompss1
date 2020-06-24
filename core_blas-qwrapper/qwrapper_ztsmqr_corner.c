/**
 *
 * @file qwrapper_ztsmqr_corner.c
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
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ztsmqr_corner(Quark *quark, Quark_Task_Flags *task_flags,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              PLASMA_Complex64_t *A1, int lda1,
                              PLASMA_Complex64_t *A2, int lda2,
                              PLASMA_Complex64_t *A3, int lda3,
                              const PLASMA_Complex64_t *V, int ldv,
                              const PLASMA_Complex64_t *T, int ldt)
{
    int ldwork = nb;

    DAG_CORE_TSMQR;
    QUARK_Insert_Task(quark, CORE_ztsmqr_corner_quark, task_flags,
        sizeof(int),                        &m1,    VALUE,
        sizeof(int),                        &n1,    VALUE,
        sizeof(int),                        &m2,    VALUE,
        sizeof(int),                        &n2,    VALUE,
        sizeof(int),                        &m3,    VALUE,
        sizeof(int),                        &n3,    VALUE,
        sizeof(int),                        &k,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        sizeof(int),                        &nb,    VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A1,            INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &lda1,  VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A2,            INOUT,
        sizeof(int),                        &lda2,  VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A3,            INOUT|QUARK_REGION_D|QUARK_REGION_L,
        sizeof(int),                        &lda3,  VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    V,             INPUT,
        sizeof(int),                        &ldv,   VALUE,
        sizeof(PLASMA_Complex64_t)*ib*nb,    T,             INPUT,
        sizeof(int),                        &ldt,   VALUE,
        sizeof(PLASMA_Complex64_t)*4*nb*nb,    NULL,          SCRATCH,
        sizeof(int),                        &ldwork, VALUE,
        0);
}


#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztsmqr_corner_quark = PCORE_ztsmqr_corner_quark
#define CORE_ztsmqr_corner_quark PCORE_ztsmqr_corner_quark
#endif
void CORE_ztsmqr_corner_quark(Quark *quark)
{
    int m1;
    int n1;
    int m2;
    int n2;
    int m3;
    int n3;
    int k;
    int ib;
    int nb;
    PLASMA_Complex64_t *A1;
    int lda1;
    PLASMA_Complex64_t *A2;
    int lda2;
    PLASMA_Complex64_t *A3;
    int lda3;
    PLASMA_Complex64_t *V;
    int ldv;
    PLASMA_Complex64_t *T;
    int ldt;
    PLASMA_Complex64_t *WORK;
    int ldwork;

    quark_unpack_args_21(quark, m1, n1, m2, n2, m3, n3, k, ib, nb,
                         A1, lda1, A2, lda2, A3, lda3, V, ldv, T, ldt, WORK, ldwork);
    CORE_ztsmqr_corner(m1, n1, m2, n2, m3, n3, k, ib, nb,
                       A1, lda1, A2, lda2, A3, lda3, V, ldv, T, ldt, WORK, ldwork);
}
