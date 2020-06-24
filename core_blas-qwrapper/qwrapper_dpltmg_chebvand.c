/**
 *
 * @file qwrapper_dpltmg_chebvand.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:44:58 2014
 *
 **/
#include "common.h"

void
CORE_dpltmg_chebvand_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dpltmg_chebvand( Quark *quark, Quark_Task_Flags *task_flags,
                                 int M, int N, double *A, int LDA,
                                 int gN, int m0, int n0,
                                 double *W )
{
    DAG_CORE_PLRNT;
    QUARK_Insert_Task(quark, CORE_dpltmg_chebvand_quark, task_flags,
        sizeof(int),                      &M,       VALUE,
        sizeof(int),                      &N,       VALUE,
        sizeof(double)*LDA*N, A,            OUTPUT,
        sizeof(int),                      &LDA,     VALUE,
        sizeof(int),                      &gN,      VALUE,
        sizeof(int),                      &m0,      VALUE,
        sizeof(int),                      &n0,      VALUE,
        sizeof(double)*N*2,   W,            INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dpltmg_chebvand_quark = PCORE_dpltmg_chebvand_quark
#define CORE_dpltmg_chebvand_quark PCORE_dpltmg_chebvand_quark
#endif
void CORE_dpltmg_chebvand_quark(Quark *quark)
{
    int m;
    int n;
    double *A;
    int lda;
    int gN;
    int m0;
    int n0;
    double *W;

    quark_unpack_args_8( quark, m, n, A, lda, gN, m0, n0, W );
    CORE_dpltmg_chebvand( m, n, A, lda, gN, m0, n0, W );
}
