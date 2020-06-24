/**
 *
 * @file core_zgeam.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeam = PCORE_zgeam
#define CORE_zgeam PCORE_zgeam
#endif
void CORE_zgeam(PLASMA_enum transA, int transB,
                int m, int n, 
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                PLASMA_Complex64_t beta, PLASMA_Complex64_t *B, int LDB)
{   
    int i, j; 
    if (transA == PlasmaNoTrans ){
        if (transB == PlasmaNoTrans){
           for (i = 0; i < m; i++){
              for (j = 0 ; j < n; j++){
                    B[j*LDB+i] = creal(alpha)*creal(A[j*LDA+i]) - cimag(alpha)*cimag(A[j*LDA+i]) +  (cimag(alpha)*creal(A[j*LDA+i]) + creal(alpha)*cimag(A[j*LDA+i]))*I +
                                 creal(beta)*creal(B[j*LDB+i]) - cimag(beta)*cimag(B[j*LDB+i]) +  (cimag(beta)*creal(B[j*LDB+i]) + creal(beta)*cimag(B[j*LDB+i]))*I;

              }
           }

        }
        else {
           for (i = 0; i < m; i++){
               for (j = 0 ; j < n; j++){
                    B[j+LDB*i] = creal(alpha)*creal(A[j*LDA+i]) - cimag(alpha)*cimag(A[j*LDA+i]) +  (cimag(alpha)*creal(A[j*LDA+i]) + creal(alpha)*cimag(A[j*LDA+i]))*I +
                                 creal(beta)*creal(B[j+LDB*i]) - cimag(beta)*cimag(B[j+LDB*i]) +  (cimag(beta)*creal(B[j+LDB*i]) + creal(beta)*cimag(B[j+LDB*i]))*I;
 
               }
           }

        }
    }
    else {
        if (transB == PlasmaNoTrans) {
           for (i = 0; i < m; i++){
               for (j = 0 ; j < n; j++){
                    B[j*LDB+i] = creal(alpha)*creal(A[j+LDA*i]) - cimag(alpha)*cimag(A[j+LDA*i]) +  (cimag(alpha)*creal(A[j+LDA*i]) + creal(alpha)*cimag(A[j+LDA*i]))*I +
                                 creal(beta)*creal(B[j*LDB+i]) - cimag(beta)*cimag(B[j*LDB+i]) +  (cimag(beta)*creal(B[j*LDB+i]) + creal(beta)*cimag(B[j*LDB+i]))*I;
  
               }
           }
       }
       else {
           for (i = 0; i < m; i++){
               for (j = 0 ; j < n; j++){
                    B[j+LDB*i] = creal(alpha)*creal(A[j+LDA*i]) - cimag(alpha)*cimag(A[j+LDA*i]) +  (cimag(alpha)*creal(A[j+LDA*i]) + creal(alpha)*cimag(A[j+LDA*i]))*I +
                                 creal(beta)*creal(B[j+LDB*i]) - cimag(beta)*cimag(B[j+LDB*i]) +  (cimag(beta)*creal(B[j+LDB*i]) + creal(beta)*cimag(B[j+LDB*i]))*I;
 
               }
           }
      }
    }
}

/***************************************************************************//**
 *
 **/

void QUARK_CORE_zgeam(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum transA, PLASMA_enum transB,
                      int m, int n, int k, int nb,
                      PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                      PLASMA_Complex64_t beta, PLASMA_Complex64_t *B, int LDB)
{
    DAG_CORE_GEAM;
    QUARK_Insert_Task(quark, CORE_zgeam_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &LDA,       VALUE,
        sizeof(PLASMA_Complex64_t),         &beta,      VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,                 INOUT,
        sizeof(int),                        &LDB,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeam_quark = PCORE_zgeam_quark
#define CORE_zgeam_quark PCORE_zgeam_quark
#endif
void CORE_zgeam_quark(Quark *quark)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    int m;
    int n;
    int k;
    int nb;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t beta;
    PLASMA_Complex64_t *B;
    int LDB;
    int i, j;
    quark_unpack_args_11(quark, transA, transB, m, n, nb,  alpha, A, LDA, beta, B, LDB);
    if (transA == PlasmaNoTrans ){
        if (transB == PlasmaNoTrans){
           for (i = 0; i < m; i++){
              for (j = 0 ; j < n; j++){
                    B[j*LDB+i] = creal(alpha)*creal(A[j*LDA+i]) - cimag(alpha)*cimag(A[j*LDA+i]) +  (cimag(alpha)*creal(A[j*LDA+i]) + creal(alpha)*cimag(A[j*LDA+i]))*I +
                                 creal(beta)*creal(B[j*LDB+i]) - cimag(beta)*cimag(B[j*LDB+i]) +  (cimag(beta)*creal(B[j*LDB+i]) + creal(beta)*cimag(B[j*LDB+i]))*I;

              }
           }

        }
        else {
           for (i = 0; i < m; i++){
               for (j = 0 ; j < n; j++){
                    B[j+LDB*i] = creal(alpha)*creal(A[j*LDA+i]) - cimag(alpha)*cimag(A[j*LDA+i]) + (cimag(alpha)*creal(A[j*LDA+i]) + creal(alpha)*cimag(A[j*LDA+i]))*I +
                                 creal(beta)*creal(B[j+LDB*i]) + cimag(beta)*cimag(B[j+LDB*i]) + (cimag(beta)*creal(B[j+LDB*i]) - creal(beta)*cimag(B[j+LDB*i]))*I;

               }
           }

        }
    }
    else {
        if (transB == PlasmaNoTrans) {
           for (i = 0; i < m; i++){
               for (j = 0 ; j < n; j++){
                    B[j*LDB+i] = creal(alpha)*creal(A[j+LDA*i]) + cimag(alpha)*cimag(A[j+LDA*i]) +  (cimag(alpha)*creal(A[j+LDA*i]) - creal(alpha)*cimag(A[j+LDA*i]))*I +
                                 creal(beta)*creal(B[j*LDB+i]) - cimag(beta)*cimag(B[j*LDB+i]) +  (cimag(beta)*creal(B[j*LDB+i]) + creal(beta)*cimag(B[j*LDB+i]))*I;

               }
           }
       }
       else {
           for (i = 0; i < m; i++){
               for (j = 0 ; j < n; j++){
                    B[j+LDB*i] = creal(alpha)*creal(A[j+LDA*i]) + cimag(alpha)*cimag(A[j+LDA*i]) +  (cimag(alpha)*creal(A[j+LDA*i]) - creal(alpha)*cimag(A[j+LDA*i]))*I +
                                 creal(beta)*creal(B[j+LDB*i]) + cimag(beta)*cimag(B[j+LDB*i]) +  (cimag(beta)*creal(B[j+LDB*i]) - creal(beta)*cimag(B[j+LDB*i]))*I;

               }
           }
      }
    }



}
