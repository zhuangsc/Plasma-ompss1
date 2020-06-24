/**
 *
 * @file example_dtrsm.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @brief Example for solving a system of linear equations using 
 *        Cholesky factorization and PLASMA_dtrsm routine
 *
 * @version 2.6.0
 * @author Bilel Hadri
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:45:21 2014
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

int check_solution(int, int, double*, int, double*, double*, int);

int IONE=1;
int ISEED[4] = {0,0,0,1};   /* initial seed for dlarnv() */

int main ()
{

    int cores = 2;
    int N     = 10 ;
    int LDA   = 10 ;
    int NRHS  = 5 ;
    int LDB   = 10 ;
    int info_solution;

    double *A1   = (double *)malloc(LDA*N*sizeof(double));
    double *A2   = (double *)malloc(LDA*N*sizeof(double));
    double *B1   = (double *)malloc(LDB*NRHS*sizeof(double));
    double *B2   = (double *)malloc(LDB*NRHS*sizeof(double));

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)||(!B1)||(!B2)){
        printf("Out of Memory \n ");
        return EXIT_SUCCESS;
    }

    /* Plasma Initialize */
    PLASMA_Init(cores);
    printf("-- PLASMA is initialized to run on %d cores. \n",cores);

    /* Initialize A1 and A2 for Symmetric Positive Matrix */
    PLASMA_dplgsy( (double)N, N, A1, LDA, 51 );
    PLASMA_dlacpy( PlasmaUpperLower, N, N, A1, LDA, A2, LDA );

    /* Initialize B1 and B2 */
    PLASMA_dplrnt( N, NRHS, B1, LDB, 371 );
    PLASMA_dlacpy( PlasmaUpperLower, N, NRHS, B1, LDB, B2, LDB );

    /* PLASMA routines */
    PLASMA_dpotrf(PlasmaLower, N, A2, LDA);
    PLASMA_dtrsm(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 
                 N, NRHS, (double)1.0, A2, LDA, B2, LDB);
    PLASMA_dtrsm(PlasmaLeft, PlasmaLower, PlasmaTrans, PlasmaNonUnit,
                 N, NRHS, (double)1.0, A2, LDA, B2, LDB);

    /* Check the solution */
    info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB);

    if ( info_solution != 0 )
       printf("-- Error in DTRSM example ! \n");
    else
       printf("-- Run of DTRSM example successful ! \n");

    free(A1); free(A2); free(B1); free(B2);

    PLASMA_Finalize();

    return EXIT_SUCCESS;
}


/*------------------------------------------------------------------------
 *  Check the accuracy of the solution of the linear system
 */

int check_solution(int N, int NRHS, double *A1, int LDA, double *B1, double *B2, int LDB)
{
    int info_solution;
    double Rnorm, Anorm, Xnorm, Bnorm;
    double alpha, beta;
    double *work = (double *)malloc(N*sizeof(double));
    double eps;

    eps = LAPACKE_dlamch_work('e');

    alpha = 1.0;
    beta  = -1.0;

    Xnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B2, LDB, work);
    Anorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, A1, LDA, work);
    Bnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, NRHS, N, (alpha), A1, LDA, B2, LDB, (beta), B1, LDB);
    Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n",Rnorm/((Anorm*Xnorm+Bnorm)*N*eps));

    if (Rnorm/((Anorm*Xnorm+Bnorm)*N*eps) > 10.0){
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
     }
    else{
        printf("-- The solution is CORRECT ! \n");
        info_solution = 0;
    }

    free(work);

    return info_solution;
}
