/**
 *
 * @file qwrapper_dswpab.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @generated d Tue Jan  7 11:44:57 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dswpab(Quark *quark, Quark_Task_Flags *task_flags,
                       int i, int n1, int n2,
                       double *A, int szeA)
{
    DAG_CORE_SWPAB;
    QUARK_Insert_Task(
        quark, CORE_dswpab_quark, task_flags,
        sizeof(int),                           &i,   VALUE,
        sizeof(int),                           &n1,  VALUE,
        sizeof(int),                           &n2,  VALUE,
        sizeof(double)*szeA,       A,            INOUT,
        sizeof(double)*min(n1,n2), NULL,         SCRATCH,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dswpab_quark = PCORE_dswpab_quark
#define CORE_dswpab_quark PCORE_dswpab_quark
#endif
void CORE_dswpab_quark(Quark *quark)
{
    int i;
    int n1;
    int n2;
    double *A;
    double *work;

    quark_unpack_args_5(quark, i, n1, n2, A, work);
    CORE_dswpab( i, n1, n2, A, work);
}
