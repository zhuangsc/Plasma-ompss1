/**
 *
 * @generated s Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_spotrf"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_POTRF( N )
#define _FADDS FADDS_POTRF( N )

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    int uplo = PlasmaLower;

    LDA = max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, float, LDA, N );

    /* Initialiaze Data */
    PLASMA_splgsy( (float)N, N, A, LDA, 51 );

    /* Save A and b  */
    PASTE_CODE_ALLOCATE_COPY( A2, check, float, A, LDA, N    );

    /* PLASMA SPOSV */
    START_TIMING();
    PLASMA_spotrf(uplo, N, A, LDA);
    STOP_TIMING();

    /* Check the solution */
    if (check)
      {
        PASTE_CODE_ALLOCATE_MATRIX( B, check, float, LDB, NRHS );
        PLASMA_splrnt( N, NRHS, B, LDB, 5673 );
        PASTE_CODE_ALLOCATE_COPY( X,  check, float, B, LDB, NRHS );

        PLASMA_spotrs(uplo, N, NRHS, A, LDA, X, LDB);

        dparam[IPARAM_RES] = s_check_solution(N, N, NRHS, A2, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]),
                                              &(dparam[IPARAM_XNORM]));

        free(A2); free(B); free(X);
      }

    free(A);

    return 0;
}
