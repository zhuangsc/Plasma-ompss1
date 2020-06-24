/**
 *
 * @file qwrapper_cher2k.c
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
 * @generated c Tue Jan  7 11:44:57 2014
 *
 **/
#include "common.h"

#undef REAL
#define COMPLEX
#ifdef COMPLEX
/***************************************************************************//**
 *
 **/
void QUARK_CORE_cher2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       PLASMA_Complex32_t alpha, const PLASMA_Complex32_t *A, int lda,
                       const PLASMA_Complex32_t *B, int ldb,
                       float beta, PLASMA_Complex32_t *C, int ldc)
{
    DAG_CORE_HER2K;
    QUARK_Insert_Task(quark, CORE_cher2k_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex32_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(float),                     &beta,      VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_cher2k_quark = PCORE_cher2k_quark
#define CORE_cher2k_quark PCORE_cher2k_quark
#endif
void CORE_cher2k_quark(Quark *quark)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    int n;
    int k;
    PLASMA_Complex32_t alpha;
    PLASMA_Complex32_t *A;
    int lda;
    PLASMA_Complex32_t *B;
    int ldb;
    float beta;
    PLASMA_Complex32_t *C;
    int ldc;

    quark_unpack_args_12(quark, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    cblas_cher2k(CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                 n, k, CBLAS_SADDR(alpha), A, lda, B, ldb, beta, C, ldc);
}
#endif
