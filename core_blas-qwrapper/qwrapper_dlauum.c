/**
 *
 * @file qwrapper_dpotrf.c
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
 * @generated d Tue Jan  7 11:44:56 2014
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dlauum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       double *A, int lda)
{
    DAG_CORE_LAUUM;
    QUARK_Insert_Task(quark, CORE_dlauum_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(double)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dlauum_quark = PCORE_dlauum_quark
#define CORE_dlauum_quark PCORE_dlauum_quark
#endif
void CORE_dlauum_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int N;
    double *A;
    int LDA;

    quark_unpack_args_4(quark, uplo, N, A, LDA);
    LAPACKE_dlauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, A, LDA);
}
