/**
 *
 * @generated d Tue Jan  7 11:45:23 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgemm"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEMM(M, N, K)
#define _FADDS FADDS_GEMM(M, N, K)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    double alpha, beta;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    LDB = max(K, iparam[IPARAM_LDB]);
    LDC = max(M, iparam[IPARAM_LDC]);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A,      1, double, LDA, K   );
	#pragma omp register([LDA*K]A)
    PASTE_CODE_ALLOCATE_MATRIX( B,      1, double, LDB, N   );
	#pragma omp register([LDB*N]B)
    PASTE_CODE_ALLOCATE_MATRIX( C,      1, double, LDC, N   );
	#pragma omp register([LDC*N]C)
    PASTE_CODE_ALLOCATE_MATRIX( C2, check, double, LDC, N   );

	int runtime = RT_get_runtime();
	int ws = RT_get_ws();
	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_QUARK);
		RT_set_ws(1);
	}

    PLASMA_dplrnt( M, K, A, LDA,  453 );
    PLASMA_dplrnt( K, N, B, LDB, 5673 );
    PLASMA_dplrnt( M, N, C, LDC,  740 );

    LAPACKE_dlarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_dlarnv_work(1, ISEED, 1, &beta );

    if (check)
    {
        memcpy(C2, C, LDC*N*sizeof(double));
    }

	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_OMPSS);
		RT_set_ws(ws);
	}

    START_TIMING();
    PLASMA_dgemm( PlasmaNoTrans, PlasmaNoTrans, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC );
    STOP_TIMING();
    
	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_QUARK);
		RT_set_ws(1);
	}

    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = d_check_gemm( PlasmaNoTrans, PlasmaNoTrans, M, N, K,
                                           alpha, A, LDA, B, LDB, beta, C, C2, LDC,
                                           &(dparam[IPARAM_ANORM]), 
                                           &(dparam[IPARAM_BNORM]), 
                                           &(dparam[IPARAM_XNORM]));
//        free(C2);
    }

//    free( A );
//    free( B );
//    free( C );

	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_OMPSS);
		RT_set_ws(ws);
	}

    return 0;
}
