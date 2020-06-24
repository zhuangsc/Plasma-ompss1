/**
 *
 * @generated s Tue Jan  7 11:45:26 2014
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include <plasma.h>
#include <core_blas.h>
#include "auxiliary.h"

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */

int s_check_orthogonality(int M, int N, int LDQ, float *Q)
{
    float alpha, beta;
    float normQ;
    int info_ortho;
    int i;
    int minMN = min(M, N);
    float eps;
    float *work = (float *)malloc(minMN*sizeof(float));

    eps = LAPACKE_slamch_work('e');
    alpha = 1.0;
    beta  = -1.0;

    /* Build the idendity matrix USE DLASET?*/
    float *Id = (float *) malloc(minMN*minMN*sizeof(float));
    memset((void*)Id, 0, minMN*minMN*sizeof(float));
    for (i = 0; i < minMN; i++)
        Id[i*minMN+i] = (float)1.0;

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_ssyrk(CblasColMajor, CblasUpper, CblasTrans, N, M, alpha, Q, LDQ, beta, Id, N);
    else
        cblas_ssyrk(CblasColMajor, CblasUpper, CblasNoTrans, M, N, alpha, Q, LDQ, beta, Id, M);

    normQ = LAPACKE_slansy_work(LAPACK_COL_MAJOR, 'i', 'u', minMN, Id, minMN, work);

    printf("============\n");
    printf("Checking the orthogonality of Q \n");
    printf("||Id-Q'*Q||_oo / (N*eps) = %e \n",normQ/(minMN*eps));

    if ( isnan(normQ / (minMN * eps)) || (normQ / (minMN * eps) > 10.0) ) {
        printf("-- Orthogonality is suspicious ! \n");
        info_ortho=1;
    }
    else {
        printf("-- Orthogonality is CORRECT ! \n");
        info_ortho=0;
    }

    free(work); free(Id);

    return info_ortho;
}

/*------------------------------------------------------------
 *  Check the factorization QR
 */

int s_check_QRfactorization(int M, int N, float *A1, float *A2, int LDA, float *Q)
{
    float Anorm, Rnorm;
    float alpha, beta;
    int info_factorization;
    int i,j;
    float eps;

    eps = LAPACKE_slamch_work('e');

    float *Ql       = (float *)malloc(M*N*sizeof(float));
    float *Residual = (float *)malloc(M*N*sizeof(float));
    float *work              = (float *)malloc(max(M,N)*sizeof(float));

    alpha=1.0;
    beta=0.0;

    if (M >= N) {
        /* Extract the R */
        float *R = (float *)malloc(N*N*sizeof(float));
        memset((void*)R, 0, N*N*sizeof(float));
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,'u', M, N, A2, LDA, R, N);

        /* Perform Ql=Q*R */
        memset((void*)Ql, 0, M*N*sizeof(float));
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, N, (alpha), Q, LDA, R, N, (beta), Ql, M);
        free(R);
    }
    else {
        /* Extract the L */
        float *L = (float *)malloc(M*M*sizeof(float));
        memset((void*)L, 0, M*M*sizeof(float));
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,'l', M, N, A2, LDA, L, M);

    /* Perform Ql=LQ */
        memset((void*)Ql, 0, M*N*sizeof(float));
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, M, (alpha), L, M, Q, LDA, (beta), Ql, M);
        free(L);
    }

    /* Compute the Residual */
    for (i = 0; i < M; i++)
        for (j = 0 ; j < N; j++)
            Residual[j*M+i] = A1[j*LDA+i]-Ql[j*M+i];

    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', M, N, Residual, M, work);
    Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', M, N, A2, LDA, work);

    if (M >= N) {
        printf("============\n");
        printf("Checking the QR Factorization \n");
        printf("-- ||A-QR||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));
    }
    else {
        printf("============\n");
        printf("Checking the LQ Factorization \n");
        printf("-- ||A-LQ||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));
    }

    if (isnan(Rnorm / (Anorm * N *eps)) || (Rnorm / (Anorm * N * eps) > 10.0) ) {
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else {
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(work); free(Ql); free(Residual);

    return info_factorization;
}

/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */

int s_check_LLTfactorization(int N, float *A1, float *A2, int LDA, int uplo)
{
    float Anorm, Rnorm;
    float alpha;
    int info_factorization;
    int i,j;
    float eps;

    eps = LAPACKE_slamch_work('e');

    float *Residual = (float *)malloc(N*N*sizeof(float));
    float *L1       = (float *)malloc(N*N*sizeof(float));
    float *L2       = (float *)malloc(N*N*sizeof(float));
    float *work              = (float *)malloc(N*sizeof(float));

    memset((void*)L1, 0, N*N*sizeof(float));
    memset((void*)L2, 0, N*N*sizeof(float));

    alpha= 1.0;

    LAPACKE_slacpy_work(LAPACK_COL_MAJOR,' ', N, N, A1, LDA, Residual, N);

    /* Dealing with L'L or U'U  */
    if (uplo == PlasmaUpper){
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L1, N);
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L2, N);
        cblas_strmm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, N, N, (alpha), L1, N, L2, N);
    }
    else{
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L1, N);
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L2, N);
        cblas_strmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, N, N, (alpha), L1, N, L2, N);
    }

    /* Compute the Residual || A -L'L|| */
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
           Residual[j*N+i] = L2[j*N+i] - Residual[j*N+i];

    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', N, N, Residual, N, work);
    Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', N, N, A1, LDA, work);

    printf("============\n");
    printf("Checking the Cholesky Factorization \n");
    printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));

    if ( isnan(Rnorm/(Anorm*N*eps)) || (Rnorm/(Anorm*N*eps) > 10.0) ){
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else{
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(Residual); free(L1); free(L2); free(work);

    return info_factorization;
}

/*--------------------------------------------------------------
 * Check the gemm
 */
float s_check_gemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K,
                   float alpha, float *A, int LDA,
                   float *B, int LDB,
                   float beta, float *Cplasma,
                   float *Cref, int LDC,
                   float *Cinitnorm, float *Cplasmanorm, float *Clapacknorm )
{
    float beta_const = -1.0;
    float Rnorm;
    float *work = (float *)malloc(max(K,max(M, N))* sizeof(float));

    *Cinitnorm   = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', M, N, Cref,    LDC, work);
    *Cplasmanorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', M, N, Cplasma, LDC, work);

    cblas_sgemm(CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K,
                (alpha), A, LDA, B, LDB, (beta), Cref, LDC);

    *Clapacknorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', M, N, Cref, LDC, work);

    cblas_saxpy(LDC * N, (beta_const), Cplasma, 1, Cref, 1);

    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', M, N, Cref, LDC, work);

    free(work);

    return Rnorm;
}

/*--------------------------------------------------------------
 * Check the trsm
 */
float s_check_trsm(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum trans, PLASMA_enum diag,
                   int M, int NRHS, float alpha,
                   float *A, int LDA,
                   float *Bplasma, float *Bref, int LDB,
                   float *Binitnorm, float *Bplasmanorm, float *Blapacknorm )
{
    float beta_const = -1.0;
    float Rnorm;
    float *work = (float *)malloc(max(M, NRHS)* sizeof(float));
    /*float eps = LAPACKE_slamch_work('e');*/

    *Binitnorm   = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', M, NRHS, Bref,    LDB, work);
    *Bplasmanorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'm', M, NRHS, Bplasma, LDB, work);

    cblas_strsm(CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag, M, NRHS,
                (alpha), A, LDA, Bref, LDB);

    *Blapacknorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'm', M, NRHS, Bref, LDB, work);

    cblas_saxpy(LDB * NRHS, (beta_const), Bplasma, 1, Bref, 1);

    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'm', M, NRHS, Bref, LDB, work);
    Rnorm = Rnorm / *Blapacknorm; 
    /* max(M,NRHS) * eps);*/

    free(work);

    return Rnorm;
}

/*--------------------------------------------------------------
 * Check the solution
 */

float s_check_solution(int M, int N, int NRHS, float *A, int LDA,
                      float *B,  float *X, int LDB,
                      float *anorm, float *bnorm, float *xnorm )
{
/*     int info_solution; */
    float Rnorm = -1.00;
    float zone  =  1.0;
    float mzone = -1.0;
    float *work = (float *)malloc(max(M, N)* sizeof(float));

    *anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', M, N,    A, LDA, work);
    *xnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', M, NRHS, X, LDB, work);
    *bnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', N, NRHS, B, LDB, work);

    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, NRHS, N, (zone), A, LDA, X, LDB, (mzone), B, LDB);

    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, 'i', N, NRHS, B, LDB, work);

    free(work);

    return Rnorm;
}
