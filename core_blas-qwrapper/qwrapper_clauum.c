/**
 *
 * @file qwrapper_cpotrf.c
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
 * @generated c Tue Jan  7 11:44:56 2014
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clauum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       PLASMA_Complex32_t *A, int lda)
{
    DAG_CORE_LAUUM;
    QUARK_Insert_Task(quark, CORE_clauum_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clauum_quark = PCORE_clauum_quark
#define CORE_clauum_quark PCORE_clauum_quark
#endif
void CORE_clauum_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int N;
    PLASMA_Complex32_t *A;
    int LDA;

    quark_unpack_args_4(quark, uplo, N, A, LDA);
    LAPACKE_clauum_work(LAPACK_COL_MAJOR, lapack_const(uplo), N, A, LDA);
}
