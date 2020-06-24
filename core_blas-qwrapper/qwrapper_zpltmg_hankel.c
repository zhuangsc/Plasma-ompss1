/**
 *
 * @file qwrapper_zpltmg_hankel.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

void
CORE_zpltmg_hankel_quark(Quark *quark);

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zpltmg_hankel( Quark *quark, Quark_Task_Flags *task_flags,
                               PLASMA_enum uplo, int M, int N, PLASMA_Complex64_t *A, int LDA,
                               int m0, int n0, int nb,
                               const PLASMA_Complex64_t *V1,
                               const PLASMA_Complex64_t *V2)
{
    DAG_CORE_PLRNT;
    QUARK_Insert_Task(quark, CORE_zpltmg_hankel_quark, task_flags,
        sizeof(PLASMA_enum),              &uplo,    VALUE,
        sizeof(int),                      &M,       VALUE,
        sizeof(int),                      &N,       VALUE,
        sizeof(PLASMA_Complex64_t)*LDA*N, A,            OUTPUT,
        sizeof(int),                      &LDA,     VALUE,
        sizeof(int),                      &m0,      VALUE,
        sizeof(int),                      &n0,      VALUE,
        sizeof(int),                      &nb,      VALUE,
        sizeof(PLASMA_Complex64_t)*nb,    V1,           INPUT,
        sizeof(PLASMA_Complex64_t)*nb,    V2,           INPUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zpltmg_hankel_quark = PCORE_zpltmg_hankel_quark
#define CORE_zpltmg_hankel_quark PCORE_zpltmg_hankel_quark
#endif
void CORE_zpltmg_hankel_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int m;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    int m0;
    int n0;
    int nb;
    const PLASMA_Complex64_t *V1;
    const PLASMA_Complex64_t *V2;

    quark_unpack_args_10( quark, uplo, m, n, A, lda, m0, n0, nb, V1, V2 );
    CORE_zpltmg_hankel( uplo, m, n, A, lda, m0, n0, nb, V1, V2 );
}
