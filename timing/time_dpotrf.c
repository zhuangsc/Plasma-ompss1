/**
 *
 * @generated d Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dpotrf"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_POTRF( N )
#define _FADDS FADDS_POTRF( N )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    int uplo = PlasmaLower;

    LDA = max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, double, LDA, N );
	#pragma omp register([LDA*N]A)
	int runtime = RT_get_runtime();
	int ws = RT_get_ws();
	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_QUARK);
		RT_set_ws(1);
	}

    /* Initialiaze Data */
    PLASMA_dplgsy( (double)N, N, A, LDA, 51 );

    /* Save A and b  */
    PASTE_CODE_ALLOCATE_COPY( A2, check, double, A, LDA, N    );

	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_OMPSS);
		RT_set_ws(ws);
	}

    /* PLASMA DPOSV */
    START_TIMING();
    PLASMA_dpotrf(uplo, N, A, LDA);
    STOP_TIMING();

	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_QUARK);
		RT_set_ws(1);
	}

    /* Check the solution */
    if (check)
      {
        PASTE_CODE_ALLOCATE_MATRIX( B, check, double, LDB, NRHS );
        PLASMA_dplrnt( N, NRHS, B, LDB, 5673 );
        PASTE_CODE_ALLOCATE_COPY( X,  check, double, B, LDB, NRHS );

        PLASMA_dpotrs(uplo, N, NRHS, A, LDA, X, LDB);

        dparam[IPARAM_RES] = d_check_solution(N, N, NRHS, A2, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]),
                                              &(dparam[IPARAM_XNORM]));

//        free(A2); free(B); free(X);
      }

//    free(A);
	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_OMPSS);
		RT_set_ws(ws);
	}

    return 0;
}
