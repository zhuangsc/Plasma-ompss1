/**
 *
 * @file example_spotrf.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @brief Example of Cholesky factorization
 *
 * @version 2.6.0
 * @author Bilel Hadri
 * @date 2010-11-15
 * @generated s Tue Jan  7 11:45:20 2014
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

int check_factorization(int, float*, float*, int, int);

int IONE=1;
int ISEED[4] = {0,0,0,1};   /* initial seed for slarnv() */

int main ()
{

    int cores = 2;
    int N     = 10 ;
    int LDA   = 10 ;
    int info_factorization;

    float *A1   = (float *)malloc(LDA*N*sizeof(float));
    float *A2   = (float *)malloc(LDA*N*sizeof(float));

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)){
        printf("Out of Memory \n ");
        return EXIT_SUCCESS;
    }

    /* Plasma Initialize */
    PLASMA_Init(cores);
    printf("-- PLASMA is initialized to run on %d cores. \n",cores);

    /* Initialize A1 and A2 for Symmetric Positive Matrix */
    PLASMA_splgsy( (float)N, N, A1, LDA, 51 );
    PLASMA_slacpy( PlasmaUpperLower, N, N, A1, LDA, A2, LDA );

    /* Plasma routines */
    PLASMA_spotrf(PlasmaUpper, N, A2, LDA);

    /* Check the factorization */
    info_factorization = check_factorization( N, A1, A2, LDA, PlasmaUpper);

    if ( info_factorization != 0 )
       printf("-- Error in SPOTRF example ! \n");
    else
       printf("-- Run of SPOTRF example successful ! \n");

    free(A1); free(A2);

    PLASMA_Finalize();

    return EXIT_SUCCESS;
}


/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */

int check_factorization(int N, float *A1, float *A2, int LDA, int uplo)
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

    Rnorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, Residual, N, work);
    Anorm = LAPACKE_slange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, A1, LDA, work);

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
