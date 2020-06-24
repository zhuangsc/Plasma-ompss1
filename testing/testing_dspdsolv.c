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

int testing_dspdsolv(int argc, char **argv)
{

    int N     = 4000;
	int LDA   = 4000;
    int info_factorization;

    double *A   = (double *)malloc(LDA*N*sizeof(double));
	double *X	 = (double *)malloc(LDA*N*sizeof(double));
	double *B	 = (double *)malloc(LDA*N*sizeof(double));

    /* Check if unable to allocate memory */
    if ((!A)||(!X)||(!B)){
        fprintf(stderr, "Out of Memory \n ");
        return 0;
    }

	PLASMA_Set( PLASMA_RUNTIME_MODE, PLASMA_OMPSS); //Toggle runtime mdoe to OmpSs
	PLASMA_Set( PLASMA_TILE_SIZE, 384); //Set up tile size

	GENMAT_SYM_FULL(N, A); //Initialize A with a SPD matrix

	LAPACKE_dlarnv(IONE, ISEED, LDA*N, X);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, N, N, 
			1.0, A, LDA, X, LDA, 0.0, B, LDA);

    PLASMA_dpotrf(PlasmaLower, N, A, LDA);
	PLASMA_dtrsm(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, N, N, 1.0, A, LDA, B, LDA);
	PLASMA_dtrsm(PlasmaLeft, PlasmaLower, PlasmaTrans, PlasmaNonUnit, N, N, 1.0, A, LDA, B, LDA);

	char infnorm = 'F';
	double Normb = LAPACKE_dlange(LAPACK_COL_MAJOR, infnorm, N, N, B, LDA);
	double Normx = LAPACKE_dlange(LAPACK_COL_MAJOR, infnorm, N, N, X, LDA);
	printf("Residual= %e, norm(b): %e, norm(x): %e\n", fabs(Normb-Normx)/Normx, Normb, Normx);

    free(A); free(X); free(B);

    return 0;
}
