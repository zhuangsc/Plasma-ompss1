/**
 *
 * @file qwrapper_dgetrf_reclap.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Piotr Luszczek
 * @date 2009-11-15
 *
 * @generated d Tue Jan  7 11:44:58 2014
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dgetrf_reclap(Quark *quark, Quark_Task_Flags *task_flags,
                              int m, int n, int nb,
                              double *A, int lda,
                              int *IPIV,
                              PLASMA_sequence *sequence, PLASMA_request *request,
                              PLASMA_bool check_info, int iinfo,
                              int nbthread)
{
    DAG_CORE_GETRF;
    QUARK_Insert_Task(quark, CORE_dgetrf_reclap_quark, task_flags,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(double)*nb*nb,    A,                     INOUT,
        sizeof(int),                        &lda,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(PLASMA_sequence*),           &sequence,      VALUE,
        sizeof(PLASMA_request*),            &request,       VALUE,
        sizeof(PLASMA_bool),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        sizeof(int),                        &nbthread,      VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dgetrf_reclap_quark = PCORE_dgetrf_reclap_quark
#define CORE_dgetrf_reclap_quark PCORE_dgetrf_reclap_quark
#endif
void CORE_dgetrf_reclap_quark(Quark* quark)
{
    int M;
    int N;
    double *A;
    int LDA;
    int *IPIV;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_bool check_info;
    int iinfo;

    int info[3];
    int maxthreads;

    quark_unpack_args_10(quark, M, N, A, LDA, IPIV, sequence, request,
                         check_info, iinfo, maxthreads );

    info[1] = QUARK_Get_RankInTask(quark);
    info[2] = maxthreads;

    CORE_dgetrf_reclap( M, N, A, LDA, IPIV, info );
    if (info[1] == 0 && info[0] != PLASMA_SUCCESS && check_info)
        plasma_sequence_flush(quark, sequence, request, iinfo + info[0] );
}
