/**
 *
 * @file example_dpotrf.c
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
 * @generated d Tue Jan  7 11:45:20 2014
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <plasma.h>
#include <cblas.h>
#include <lapacke.h>
#include <core_blas.h>
#include "testing_dmain.h"

static int check_factorization(int, double*, double*, int, int);

//int IONE=1;
//int ISEED[4] = {0,0,0,1};   /* initial seed for dlarnv() */

static void GENMAT_SYM_FULL(int m, double *A) 
{
	srand48(time(NULL));

	int j;
	for (j = 0; j < m; ++j ) {
		int i;
		for( i = j; i < m; ++i ) {
			double dran = drand48();
			A[j*m+i] = A[i*m+j] = dran;
		}
  	}
	for(j = 0; j < m; ++j)
		A[j*m+j] += 10 * m;
}

int testing_dpotrf(int argc, char **argv)
{

	if (argc != 2){
		fprintf(stderr,"POTRF: N LDA\n");
		return -1;
	}

    int N     = atoi(argv[0]);
    int LDA   = atoi(argv[1]);
    int info_factorization;

	PLASMA_Set(PLASMA_TILE_SIZE, 388);
    double *A1   = (double *)malloc(LDA*N*sizeof(double));
	#pragma omp register ([LDA*N]A1)
    double *A2   = (double *)malloc(LDA*N*sizeof(double));
	#pragma omp register ([LDA*N]A2)

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)){
        printf("Out of Memory \n ");
        return 0;
    }

    /* Initialize A1 and A2 for Symmetric Positive Matrix */
	GENMAT_SYM_FULL(N, A1);
	int i;
	for(i = 0; i < N*LDA; ++i){
		A2[i] = A1[i];
	}

    /* Plasma routines */
    PLASMA_dpotrf(PlasmaUpper, N, A2, LDA);

    /* Check the factorization */
    info_factorization = check_factorization( N, A1, A2, LDA, PlasmaUpper);

    if ( info_factorization != 0 )
       printf("-- Error in DPOTRF example ! \n");
    else
       printf("-- Run of DPOTRF example successful ! \n");

    free(A1); free(A2);

    return 0;
}

static int check_factorization(int N, double *A1, double *A2, int LDA, int uplo)
{
    double Anorm, Rnorm;
    double alpha;
    int info_factorization;
    int i,j;
    double eps;

    eps = LAPACKE_dlamch_work('e');

    double *Residual = (double *)malloc(N*N*sizeof(double));
    double *L1       = (double *)malloc(N*N*sizeof(double));
    double *L2       = (double *)malloc(N*N*sizeof(double));
    double *work              = (double *)malloc(N*sizeof(double));

    memset((void*)L1, 0, N*N*sizeof(double));
    memset((void*)L2, 0, N*N*sizeof(double));

    alpha= 1.0;

    LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,' ', N, N, A1, LDA, Residual, N);

    /* Dealing with L'L or U'U  */
    if (uplo == PlasmaUpper){
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L1, N);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L2, N);
        cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, N, N, (alpha), L1, N, L2, N);
    }
    else{
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L1, N);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L2, N);
        cblas_dtrmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, N, N, (alpha), L1, N, L2, N);
    }

    /* Compute the Residual || A -L'L|| */
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
           Residual[j*N+i] = L2[j*N+i] - Residual[j*N+i];

    Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, Residual, N, work);
    Anorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, A1, LDA, work);

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

int timing_dpotrf(int argc, char **argv)
{

    int N;
    int LDA;
	int bs;
	int rep;

	double start, end;
	double elapsed;
	FILE *log;
	int num_threads;

	if (argc != 4 ) {
		fprintf(stderr,"DPOTRF n lda bs rep\n");
		exit(1);
	}

	sscanf(argv[0], "%d", &N);
	sscanf(argv[1], "%d", &LDA);
	sscanf(argv[2], "%d", &bs);
	sscanf(argv[3], "%d", &rep);

    double *A = malloc(LDA*N*sizeof(double));
	#pragma omp register ([LDA*N]A)

    /* Check if unable to allocate memory */
    if (!A){
        printf("Out of Memory \n ");
        return 0;
    }

    /* Plasma Initialize */
	PLASMA_Set( PLASMA_TILE_SIZE, bs);

	GENMAT_SYM_FULL(N, A);

	elapsed = 0.0;
	int i;
	for ( i = 0; i < rep; i++ ) {
		start = gtime();
		PLASMA_dpotrf(PlasmaUpper, N, A, LDA);
		end = gtime();
		elapsed += end - start;
	}

	num_threads = omp_get_max_threads();
	dump_info("plasma_dpotrf.log", num_threads, elapsed, rep);

    free(A);

    return 0;

}
