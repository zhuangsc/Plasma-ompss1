/**
 *
 * @file qwrapper_dpamm.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Dulceneia Becker
 * @date 2011-06-14
 * @generated d Tue Jan  7 11:44:59 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void
QUARK_CORE_dpamm(Quark *quark, Quark_Task_Flags *task_flags,
                 int op, PLASMA_enum side, int storev,
                 int m, int n, int k, int l,
                 const double *A1, int lda1,
                       double *A2, int lda2,
                 const double *V, int ldv,
                       double *W, int ldw)
{
    QUARK_Insert_Task(quark, CORE_dpamm_quark, task_flags,
        sizeof(int),                        &op,      VALUE,
        sizeof(PLASMA_enum),                &side,    VALUE,
        sizeof(PLASMA_enum),                &storev,  VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(int),                        &k,       VALUE,
        sizeof(int),                        &l,       VALUE,
        sizeof(double)*m*k,     A1,           INPUT,
        sizeof(int),                        &lda1,    VALUE,
        sizeof(double)*k*n,     A2,           INOUT,
        sizeof(int),                        &lda2,    VALUE,
        sizeof(double)*m*n,     V,            INPUT,
        sizeof(int),                        &ldv,     VALUE,
        sizeof(double)*m*n,     W,            INOUT,
        sizeof(int),                        &ldw,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
void
CORE_dpamm_quark(Quark *quark)
{
    int op;
    PLASMA_enum side;
    int storev;
    int M;
    int N;
    int K;
    int L;
    double *A1;
    int LDA1;
    double *A2;
    int LDA2;
    double *V;
    int LDV;
    double *W;
    int LDW;

    quark_unpack_args_15(quark, op, side, storev, M, N, K, L,
            A1, LDA1, A2, LDA2, V, LDV, W, LDW);

    CORE_dpamm( op, side, storev, M, N, K, L, A1, LDA1, A2, LDA2, V, LDV, W, LDW);
}
