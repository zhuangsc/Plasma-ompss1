/**
 *
 * @generated s Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sposv"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( N ) + FMULS_POTRS( N, NRHS ))
#define _FADDS (FADDS_POTRF( N ) + FADDS_POTRS( N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_enum uplo = PlasmaUpper;
    
    LDA = max(LDA, N);
    
    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, float, LDA, N    );
    PASTE_CODE_ALLOCATE_MATRIX( X, 1, float, LDB, NRHS );
    
    /* Initialiaze Data */
    PLASMA_splgsy((float)N, N, A, LDA, 51 );
    PLASMA_splrnt( N, NRHS, X, LDB, 5673 );

    /* Save A and b  */
    PASTE_CODE_ALLOCATE_COPY( Acpy, check, float, A, LDA, N    );
    PASTE_CODE_ALLOCATE_COPY( B,    check, float, X, LDB, NRHS );

    /* PLASMA SPOSV */
    START_TIMING();
    PLASMA_sposv(uplo, N, NRHS, A, LDA, X, LDB);
    STOP_TIMING();

    /* Check the solution */
    if (check)
      {
        dparam[IPARAM_RES] = s_check_solution(N, N, NRHS, Acpy, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(Acpy); free(B);
      }

    free(A); free(X); 

    return 0;
}
