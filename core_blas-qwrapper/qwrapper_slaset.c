/**
 *
 * @file qwrapper_slaset.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:44:58 2014
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_slaset(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int M, int N,
                       float alpha, float beta,
                       float *A, int LDA)
{
    DAG_CORE_LASET;
    QUARK_Insert_Task(quark, CORE_slaset_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(float),         &alpha, VALUE,
        sizeof(float),         &beta,  VALUE,
        sizeof(float)*LDA*N,    A,      OUTPUT,
        sizeof(int),                        &LDA,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_slaset_quark = PCORE_slaset_quark
#define CORE_slaset_quark PCORE_slaset_quark
#endif
void CORE_slaset_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    float alpha;
    float beta;
    float *A;
    int LDA;

    quark_unpack_args_7(quark, uplo, M, N, alpha, beta, A, LDA);
    LAPACKE_slaset_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        M, N, alpha, beta, A, LDA);
}
