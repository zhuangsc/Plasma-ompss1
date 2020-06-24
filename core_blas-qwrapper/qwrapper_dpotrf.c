/**
 *
 * @file qwrapper_dpotrf.c
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
 * @generated d Tue Jan  7 11:44:56 2014
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dpotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       double *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
    DAG_CORE_POTRF;
    QUARK_Insert_Task(quark, CORE_dpotrf_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(double)*nb*nb,    A,                 INOUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_sequence*),           &sequence,  VALUE,
        sizeof(PLASMA_request*),            &request,   VALUE,
        sizeof(int),                        &iinfo,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dpotrf_quark = PCORE_dpotrf_quark
#define CORE_dpotrf_quark PCORE_dpotrf_quark
#endif
void CORE_dpotrf_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int n;
    double *A;
    int lda;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;

    int info;

    quark_unpack_args_7(quark, uplo, n, A, lda, sequence, request, iinfo);
    info = LAPACKE_dpotrf_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        n, A, lda);
    if (sequence->status == PLASMA_SUCCESS && info != 0)
        plasma_sequence_flush(quark, sequence, request, iinfo+info);
}
