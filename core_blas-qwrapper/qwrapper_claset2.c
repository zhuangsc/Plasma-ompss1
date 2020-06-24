/**
 *
 * @file qwrapper_claset2.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated c Tue Jan  7 11:44:59 2014
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_claset2(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int M, int N,
                       PLASMA_Complex32_t alpha, PLASMA_Complex32_t *A, int LDA)
{
    DAG_CORE_LASET;
    QUARK_Insert_Task(quark, CORE_claset2_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex32_t)*M*N,     A,      OUTPUT,
        sizeof(int),                        &LDA,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_claset2_quark = PCORE_claset2_quark
#define CORE_claset2_quark PCORE_claset2_quark
#endif
void CORE_claset2_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int LDA;

    quark_unpack_args_6(quark, uplo, M, N, alpha, A, LDA);
    CORE_claset2(uplo, M, N, alpha, A, LDA);
}
