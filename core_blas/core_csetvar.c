/**
 *
 * @file core_csetvar.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mark Gates
 * @date 2010-11-15
 * @generated c Tue Jan  7 11:44:49 2014
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex32_t
 *
 * CORE_csetvar sets a single variable, x := alpha.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *         Scalar to set x to, passed by pointer so it can depend on runtime value.
 *
 * @param[out] x
 *         On exit, x = alpha.
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_csetvar = PCORE_csetvar
#define CORE_csetvar PCORE_csetvar
#endif
void CORE_csetvar(const PLASMA_Complex32_t *alpha, PLASMA_Complex32_t *x)
{
    *x = *alpha;
}
