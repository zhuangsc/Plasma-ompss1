/**
 *
 * @file qwrapper_zlacpy_pivot.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2013-02-01
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlacpy_pivot(Quark *quark, Quark_Task_Flags *task_flags,
                             const PLASMA_desc descA,
                             PLASMA_enum direct, int k1, int k2, const int *ipiv,
                             int *rankin, int *rankout,
                             PLASMA_Complex64_t *A, int lda,
                             int pos, int init)
{
    DAG_SET_PROPERTIES( "CPY_PIV"  , "white"   );
    QUARK_Insert_Task(quark, CORE_zlacpy_pivot_quark, task_flags,
        sizeof(PLASMA_desc),                    &descA,         VALUE,
        sizeof(PLASMA_enum),                    &direct,        VALUE,
        sizeof(int),                            &k1,            VALUE,
        sizeof(int),                            &k2,            VALUE,
        sizeof(int)*lda,                         ipiv,                INPUT,
        sizeof(int)*lda,                         rankin,              INOUT,
        sizeof(int)*lda,                         rankout,             OUTPUT | GATHERV,
        sizeof(PLASMA_Complex64_t)*lda*descA.nb, A,                   INOUT | GATHERV,
        sizeof(int),                            &lda,           VALUE,
        sizeof(int),                            &pos,           VALUE,
        sizeof(int),                            &init,          VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlacpy_pivot_quark = PCORE_zlacpy_pivot_quark
#define CORE_zlacpy_pivot_quark PCORE_zlacpy_pivot_quark
#endif
void CORE_zlacpy_pivot_quark(Quark *quark)
{
    PLASMA_desc descA;
    PLASMA_enum direct;
    PLASMA_Complex64_t *A;
    int lda, pos, k1, k2;
    int *rankin;
    int *rankout;
    const int *ipiv;
    int init;

    quark_unpack_args_11(quark, descA, direct, k1, k2, ipiv, rankin, rankout, A, lda, pos, init);
    A = (direct == PlasmaRowwise) ? A + pos : A + (pos * lda);
    CORE_zlacpy_pivot(descA, direct, k1, k2, ipiv, rankin, rankout+pos, A, lda, init );
}
