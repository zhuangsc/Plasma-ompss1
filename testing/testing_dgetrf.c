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

static int check_factorization(int, int, double*, double*, int, int*);

int testing_dgetrf(int argc, char **argv)
{
	int M	  = 1000;
    int N     = 500;
    int LDA   = 1000;
    int info_factorization;

    double *A1   = (double *)malloc(LDA*N*sizeof(double));
    double *A2   = (double *)malloc(LDA*N*sizeof(double));

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)){
        printf("Out of Memory \n ");
        return 0;
    }

    /* Initialize A1 and A2 */
	LAPACKE_dlarnv(3,ISEED, LDA*N, A1);

	int i;
	for(i = 0; i < N*LDA; ++i){
		A2[i] = A1[i];
	}
    /* Plasma routines */
	int *IPIV = calloc(LDA, sizeof(int));
    PLASMA_dgetrf(M, N, A2, LDA, IPIV);

    /* Check the factorization */
    info_factorization = check_factorization( M, N, A1, A2, LDA, IPIV);

    if ( info_factorization != 0 )
       printf("-- Error in DGETRF example ! \n");
    else
       printf("-- Run of DGETRF example successful ! \n");

    free(A1); free(A2);

    return 0;
}

static int check_factorization(int M, int N, double *A1, double *A2, int LDA, int *IPIV)
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
	int *ipiv0 = calloc(LDA, sizeof(int));
	int dff = 0;

	LAPACKE_dgetrf(LAPACK_COL_MAJOR, M, N, A1, M, ipiv0);

	for(j=0; j<N; j++)
		dff += IPIV[j] - ipiv0[j];
	printf("dff: %d\n", dff);

    Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, A1, M, work);
    Anorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, A2, LDA, work);

	Rnorm1 = PLASMA_dlange(PlasmaInfNorm, M, N, A1, LDA);
	Anorm1 = PLASMA_dlange(PlasmaInfNorm, M, N, A2, LDA);

	printf("|Rnorm-Rnorm1|: %e, |Anorm-Anorm1|: %e\n", fabs(Rnorm-Rnorm1), fabs(Anorm-Anorm1));
    printf("============\n");
    printf("Checking the LU Factorization \n");
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
