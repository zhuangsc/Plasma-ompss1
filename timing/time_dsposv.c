/**
 *
 * @generated ds Tue Jan  7 11:45:23 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dposv"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( N ) + FMULS_POTRS( N, NRHS ))
#define _FADDS (FADDS_POTRF( N ) + FADDS_POTRS( N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    int iter;

    LDA = max(LDA, N);
    
    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, double, LDA, N    );
    PASTE_CODE_ALLOCATE_MATRIX( B, 1, double, LDB, NRHS );

    /* Initialiaze Data */
    PLASMA_dplgsy((double)N, N, A, LDA, 51 );
    PLASMA_dplrnt( N, NRHS, B, LDB, 5673 );

    PASTE_CODE_ALLOCATE_COPY( X, 1, double, B, LDB, NRHS );

    /* PLASMA DSPOSV */
    START_TIMING();
    PLASMA_dsposv(PlasmaUpper, N, NRHS, A, LDA, B, LDB, X, LDB, &iter);
    STOP_TIMING();

    /* Check the solution */
    if (check)
      {
        dparam[IPARAM_RES] = d_check_solution(N, N, NRHS, A, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
      }

    free(A); free(B); free(X); 
    return 0;
}
