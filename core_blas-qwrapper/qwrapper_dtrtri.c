/**
 *
 * @file qwrapper_dtrtri.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:44:56 2014
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dtrtri(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum diag,
                       int n, int nb,
                       double *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
    QUARK_Insert_Task(
        quark, CORE_dtrtri_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
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
#pragma weak CORE_dtrtri_quark = PCORE_dtrtri_quark
#define CORE_dtrtri_quark PCORE_dtrtri_quark
#endif
void CORE_dtrtri_quark(Quark *quark)
{
    PLASMA_enum uplo;
    PLASMA_enum diag;
    int N;
    double *A;
    int LDA;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;

    int info;

    quark_unpack_args_8(quark, uplo, diag, N, A, LDA, sequence, request, iinfo);
    info = LAPACKE_dtrtri_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo), lapack_const(diag),
        N, A, LDA);
    if ((sequence->status == PLASMA_SUCCESS) && (info > 0))
        plasma_sequence_flush(quark, sequence, request, iinfo + info);
}
