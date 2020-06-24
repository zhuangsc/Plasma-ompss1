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

static int check_factorization(int, int, double*, double*, int);

static void GENMAT_SYM_FULL(int m, int n, double *A) 
{
	srand48(time(NULL));

	int j;
	for (j = 0; j < m; ++j ) {
		int i;
		for( i = j; i < n; ++i ) {
			double dran = drand48();
			A[j*m+i] = A[i*m+j] = dran;
		}
  	}
	for(j = 0; j < m; ++j)
		A[j*m+j] += 10 * m;
}

int testing_dgeqrf(int argc, char **argv)
{

	int M	  = 1000;
    int N     = 1000;
    int LDA   = 1000;
    int info_factorization;

    double *A1   = (double *)malloc(LDA*N*sizeof(double));
    double *A2   = (double *)malloc(LDA*N*sizeof(double));
	#pragma omp register ([LDA*N]A2)

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)){
        printf("Out of Memory \n ");
        return EXIT_SUCCESS;
    }

    /* Initialize A1 and A2 for Symmetric Positive Matrix */
	GENMAT_SYM_FULL(LDA, N, A1);

	int i;
	for(i = 0; i < N*LDA; ++i){
		A2[i] = A1[i];
	}
    /* Plasma routines */
	PLASMA_desc *T;
	PLASMA_Alloc_Workspace_dgels(M, N, &T);
    PLASMA_dgeqrf(M, N, A2, LDA, T);

    /* Check the factorization */
    info_factorization = check_factorization( M, N, A1, A2, LDA);

    if ( info_factorization != 0 )
       printf("-- Error in DPOTRF example ! \n");
    else
       printf("-- Run of DPOTRF example successful ! \n");

    free(A1); free(A2);

    return 0;
}

static int check_factorization(int M, int N, double *A1, double *A2, int LDA)
{
    double Anorm, Rnorm;
	double Anorm1, Rnorm1;
    double alpha;
    int info_factorization;
    int i,j;
    double eps;

    eps = LAPACKE_dlamch_work('e');

//    double *Residual = (double *)malloc(M*N*sizeof(double));
//    double *L1       = (double *)malloc(M*N*sizeof(double));
//    double *L2       = (double *)malloc(M*N*sizeof(double));
    double *work              = (double *)malloc(N*sizeof(double));

	LAPACKE_dgeqrf(LAPACK_COL_MAJOR, M, N, A1, M, work);

    Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, A1, M, work);
    Anorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, A2, LDA, work);

	Rnorm1 = PLASMA_dlange(PlasmaInfNorm, M, N, A1, LDA);
	Anorm1 = PLASMA_dlange(PlasmaInfNorm, M, N, A2, LDA);

	printf("|Rnorm-Rnorm1|: %e, |Anorm-Anorm1|: %e\n", fabs(Rnorm-Rnorm1), fabs(Anorm-Anorm1));
    printf("============\n");
    printf("Checking the QR Factorization \n");
    printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n",fabs(Rnorm-Anorm)/(Anorm));

    if ( isnan(fabs(Rnorm-Anorm)/(Anorm)) || (fabs(Rnorm-Anorm)/(Anorm) > 10.0) ){
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else{
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

//    free(Residual); free(L1); free(L2); free(work);

    return info_factorization;
}
