/**
 *
 * @file qwrapper_csetvar.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated c Tue Jan  7 11:44:59 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * QUARK_CORE_csetvar sets a single variable, x := alpha.
 * Since x can be in the middle of a tile, we need to depend on the whole tile,
 * so add xlock argument.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *         Scalar to set x to, passed by pointer so it can depend on runtime value.
 *
 * @param[out] x
 *         On exit, x = alpha.
 *
 * @param[out] xlock
 *         Pointer to tile owning output variable x.
 *
 **/
void QUARK_CORE_csetvar(Quark *quark, Quark_Task_Flags *task_flags,
                        const PLASMA_Complex32_t *alpha, PLASMA_Complex32_t *x,
                        PLASMA_Complex32_t *xlock)
{
    DAG_SET_PROPERTIES("setvar", "orange");
    QUARK_Insert_Task(quark, CORE_csetvar_quark, task_flags,
        sizeof(PLASMA_Complex32_t),  alpha,  INPUT,
        sizeof(PLASMA_Complex32_t),  x,              NODEP,  /* INOUT; see xlock */
        sizeof(PLASMA_Complex32_t),  xlock,          INOUT,  /* tile dependency containing x */
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_csetvar_quark = PCORE_csetvar_quark
#define CORE_csetvar_quark PCORE_csetvar_quark
#endif
void CORE_csetvar_quark(Quark *quark)
{
    const PLASMA_Complex32_t *alpha;
    PLASMA_Complex32_t *x;

    quark_unpack_args_2( quark, alpha, x );
    *x = *alpha;
}
