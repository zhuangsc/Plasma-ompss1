/**
 *
 * @file qwrapper_strsm.c
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
 * @generated s Tue Jan  7 11:44:56 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_strsm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      float alpha, const float *A, int lda,
                      float *B, int ldb)
{
    DAG_CORE_TRSM;
    QUARK_Insert_Task(quark, CORE_strsm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(float),         &alpha,     VALUE,
        sizeof(float)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(float)*nb*nb,    B,                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_strsm_quark = PCORE_strsm_quark
#define CORE_strsm_quark PCORE_strsm_quark
#endif
void CORE_strsm_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum uplo;
    PLASMA_enum transA;
    PLASMA_enum diag;
    int m;
    int n;
    float alpha;
    float *A;
    int lda;
    float *B;
    int ldb;

    quark_unpack_args_11(quark, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
    cblas_strsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        m, n,
        (alpha), A, lda,
        B, ldb);
}
