/**
 *
 * @generated d Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgetrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(M, N)
#define _FADDS FADDS_GETRF(M, N)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

	int runtime = RT_get_runtime();
	int ws = RT_get_ws();
	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_QUARK);
		RT_set_ws(1);
	}
    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, double, PlasmaRealDouble, LDA, M, N );
	double *mat_ptr = descA->mat;
	#pragma omp register ( [LDA*N]mat_ptr)
//	printf("register: mat: %p, size: %d\n", mat_ptr, LDA*N);
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, min(M, N), 1 );

//	RT_runtime_info();
    PLASMA_dplrnt_Tile(descA, 3456);

    /* Save AT in lapack layout for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, double, LDA, N );

	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_OMPSS);
		RT_set_ws(ws);
	}

//	RT_runtime_info();
    START_TIMING();
    PLASMA_dgetrf_Tile( descA, piv );
    STOP_TIMING();

	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_QUARK);
		RT_set_ws(1);
	}

    /* Check the solution */
    if ( check )
    {
        PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, double, PlasmaRealDouble, LDB, N, NRHS );
        PLASMA_dplrnt_Tile( descB, 7732 );
        PASTE_TILE_TO_LAPACK( descB, b, check, double, LDB, NRHS );

        PLASMA_dgetrs_Tile( PlasmaNoTrans, descA, piv, descB );

        PASTE_TILE_TO_LAPACK( descB, x, check, double, LDB, NRHS );
        dparam[IPARAM_RES] = d_check_solution(M, N, NRHS, A, LDA, b, x, LDB,
                                              &(dparam[IPARAM_ANORM]),
                                              &(dparam[IPARAM_BNORM]),
                                              &(dparam[IPARAM_XNORM]));

//        PASTE_CODE_FREE_MATRIX( descB );
//        free(A); free(b); free(x);
    }

//	RT_set_ws(ps);
//    PASTE_CODE_FREE_MATRIX( descA );
//    free( piv );
	if ( runtime == PLASMA_OMPSS ) {
		PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_OMPSS);
		RT_set_ws(ws);
	}

    return 0;
}
