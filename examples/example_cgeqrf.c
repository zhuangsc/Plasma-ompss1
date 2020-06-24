/**
 *
 * @file example_cgeqrf.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @brief Example using QR factorization
 *
 * @version 2.6.0
 * @author Bilel Hadri
 * @date 2010-11-15
 * @generated c Tue Jan  7 11:45:21 2014
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <plasma.h>
#include <cblas.h>
#include <lapacke.h>
#include <core_blas.h>

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

int check_orthogonality(int, int, int, PLASMA_Complex32_t*);
int check_factorization(int, int, PLASMA_Complex32_t*, PLASMA_Complex32_t*, int, PLASMA_Complex32_t*);

int IONE=1;
int ISEED[4] = {0,0,0,1};   /* initial seed for clarnv() */

int main ()
{

    int cores = 2;
    int M     = 15;
    int N     = 10;
    int LDA   = 15;
    int K = min(M, N);
    int info;
    int info_ortho, info_factorization;
    int i,j;
    int LDAxN = LDA*N;

    PLASMA_Complex32_t *A1 = (PLASMA_Complex32_t *)malloc(LDA*N*sizeof(PLASMA_Complex32_t));
    PLASMA_Complex32_t *A2 = (PLASMA_Complex32_t *)malloc(LDA*N*sizeof(PLASMA_Complex32_t));
    PLASMA_Complex32_t *Q  = (PLASMA_Complex32_t *)malloc(LDA*N*sizeof(PLASMA_Complex32_t));
    PLASMA_desc *T;

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)||(!Q)){
        printf("Out of Memory \n ");
        return EXIT_SUCCESS;
    }

    /* Plasma Initialization */
    PLASMA_Init(cores);
    printf("-- PLASMA is initialized to run on %d cores. \n",cores);

    /* Allocate T */
    PLASMA_Alloc_Workspace_cgeqrf(M, N, &T);

    /* Initialize A1 and A2 */
    LAPACKE_clarnv_work(IONE, ISEED, LDAxN, A1);
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i] ;

    /* Factorization QR of the matrix A2 */
    info = PLASMA_cgeqrf(M, N, A2, LDA, T);

    /* Building the economy-size Q */
    memset((void*)Q, 0, LDA*N*sizeof(PLASMA_Complex32_t));
    for (i = 0; i < K; i++)
        Q[LDA*i+i] = 1.0;

    PLASMA_cungqr(M, N, K, A2, LDA, T, Q, LDA);

    /* Check the orthogonality, factorization and the solution */
    info_ortho = check_orthogonality(M, N, LDA, Q);
    info_factorization = check_factorization(M, N, A1, A2, LDA, Q);
    printf("--- info %d %d %d \n",info_factorization,info_ortho,info);
    if ((info_ortho != 0)|(info_factorization != 0)|(info != 0))
       printf("-- Error in CGEQRF example ! \n");
    else
       printf("-- Run of CGEQRF example successful ! \n");

    free(A1); free(A2); free(Q); free(T);

    PLASMA_Finalize();

    return EXIT_SUCCESS;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */

int check_orthogonality(int M, int N, int LDQ, PLASMA_Complex32_t *Q)
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
    PLASMA_Complex32_t *Id = (PLASMA_Complex32_t *) malloc(minMN*minMN*sizeof(PLASMA_Complex32_t));
    memset((void*)Id, 0, minMN*minMN*sizeof(PLASMA_Complex32_t));
    for (i = 0; i < minMN; i++)
        Id[i*minMN+i] = (PLASMA_Complex32_t)1.0;

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_cherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, alpha, Q, LDQ, beta, Id, N);
    else
        cblas_cherk(CblasColMajor, CblasUpper, CblasNoTrans, M, N, alpha, Q, LDQ, beta, Id, M);

    normQ = LAPACKE_clansy_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), 'u', minMN, Id, minMN, work);

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

int check_factorization(int M, int N, PLASMA_Complex32_t *A1, PLASMA_Complex32_t *A2, int LDA, PLASMA_Complex32_t *Q)
{
    float Anorm, Rnorm;
    PLASMA_Complex32_t alpha, beta;
    int info_factorization;
    int i,j;
    float eps;

    eps = LAPACKE_slamch_work('e');

    PLASMA_Complex32_t *Ql       = (PLASMA_Complex32_t *)malloc(M*N*sizeof(PLASMA_Complex32_t));
    PLASMA_Complex32_t *Residual = (PLASMA_Complex32_t *)malloc(M*N*sizeof(PLASMA_Complex32_t));
    float *work              = (float *)malloc(max(M,N)*sizeof(float));

    alpha=1.0;
    beta=0.0;

    if (M >= N) {
        /* Extract the R */
        PLASMA_Complex32_t *R = (PLASMA_Complex32_t *)malloc(N*N*sizeof(PLASMA_Complex32_t));
        memset((void*)R, 0, N*N*sizeof(PLASMA_Complex32_t));
        LAPACKE_clacpy_work(LAPACK_COL_MAJOR,'u', M, N, A2, LDA, R, N);

        /* Perform Ql=Q*R */
        memset((void*)Ql, 0, M*N*sizeof(PLASMA_Complex32_t));
        cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, N, CBLAS_SADDR(alpha), Q, LDA, R, N, CBLAS_SADDR(beta), Ql, M);
        free(R);
    }
    else {
        /* Extract the L */
        PLASMA_Complex32_t *L = (PLASMA_Complex32_t *)malloc(M*M*sizeof(PLASMA_Complex32_t));
        memset((void*)L, 0, M*M*sizeof(PLASMA_Complex32_t));
        LAPACKE_clacpy_work(LAPACK_COL_MAJOR,'l', M, N, A2, LDA, L, M);

    /* Perform Ql=LQ */
        memset((void*)Ql, 0, M*N*sizeof(PLASMA_Complex32_t));
        cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, M, CBLAS_SADDR(alpha), L, M, Q, LDA, CBLAS_SADDR(beta), Ql, M);
        free(L);
    }

    /* Compute the Residual */
    for (i = 0; i < M; i++)
        for (j = 0 ; j < N; j++)
            Residual[j*M+i] = A1[j*LDA+i]-Ql[j*M+i];

    Rnorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, Residual, M, work);
    Anorm = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, A2, LDA, work);

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
