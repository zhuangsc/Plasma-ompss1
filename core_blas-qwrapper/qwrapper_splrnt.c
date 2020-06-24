/**
 *
 * @file qwrapper_splrnt.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:44:58 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_splrnt( Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, float *A, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed )
{
    DAG_CORE_PLRNT;
    QUARK_Insert_Task(quark, CORE_splrnt_quark, task_flags,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(float)*lda*n, A,         OUTPUT,
        sizeof(int),                      &lda,  VALUE,
        sizeof(int),                      &bigM, VALUE,
        sizeof(int),                      &m0,   VALUE,
        sizeof(int),                      &n0,   VALUE,
        sizeof(unsigned long long int),   &seed, VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_splrnt_quark = PCORE_splrnt_quark
#define CORE_splrnt_quark PCORE_splrnt_quark
#endif
void CORE_splrnt_quark(Quark *quark)
{
    int m;
    int n;
    float *A;
    int lda;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    quark_unpack_args_8( quark, m, n, A, lda, bigM, m0, n0, seed );
    CORE_splrnt( m, n, A, lda, bigM, m0, n0, seed );
}

