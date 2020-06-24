/**
 *
 * @generated d Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgemm_Tile"
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
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, double, PlasmaRealDouble, LDA, M, K );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, double, PlasmaRealDouble, LDB, K, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descC, 1, double, PlasmaRealDouble, LDC, M, N );
	double *ptra = descA->mat;
	#pragma omp register([LDA*K]ptra)
	double *ptrb = descB->mat;
	#pragma omp register([LDB*N]ptrb)
	double *ptrc = descC->mat;
	#pragma omp register([LDC*N]ptrc)

	int runtime = RT_get_runtime();
	int ws = RT_get_ws();
	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_QUARK);
		RT_set_ws(1);
	}

    /* Initialiaze Data */
    PLASMA_dplrnt_Tile( descA, 5373 );
    PLASMA_dplrnt_Tile( descB, 7672 );
    PLASMA_dplrnt_Tile( descC, 6387 );
    
    LAPACKE_dlarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_dlarnv_work(1, ISEED, 1, &beta);
    
    /* Save C for check */
    PASTE_TILE_TO_LAPACK( descC, C2, check, double, LDC, N );

	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_OMPSS);
		RT_set_ws(ws);
	}

    START_TIMING();
    PLASMA_dgemm_Tile( PlasmaNoTrans, PlasmaNoTrans, alpha, descA, descB, beta, descC );
    STOP_TIMING();
    
	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_QUARK);
		RT_set_ws(1);
	}

    /* Check the solution */
    if (check)
    {
        PASTE_TILE_TO_LAPACK( descA, A, check, double, LDA, K );
        PASTE_TILE_TO_LAPACK( descB, B, check, double, LDB, N );
        PASTE_TILE_TO_LAPACK( descC, C, check, double, LDC, N );

        dparam[IPARAM_RES] = d_check_gemm( PlasmaNoTrans, PlasmaNoTrans, M, N, K,
                                           alpha, A, LDA, B, LDB, beta, C, C2, LDC,
                                           &(dparam[IPARAM_ANORM]), 
                                           &(dparam[IPARAM_BNORM]), 
                                           &(dparam[IPARAM_XNORM]));
//        free(A); free(B); free(C); free(C2);
    }

//    PASTE_CODE_FREE_MATRIX( descA );
//    PASTE_CODE_FREE_MATRIX( descB );
//    PASTE_CODE_FREE_MATRIX( descC );
	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_OMPSS);
		RT_set_ws(ws);
	}
    return 0;
}
