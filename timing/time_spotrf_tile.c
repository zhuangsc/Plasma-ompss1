/**
 *
 * @generated s Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_spotrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_POTRF( N )
#define _FADDS FADDS_POTRF( N )

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    int uplo = PlasmaUpper;

    LDA = max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1,     float, PlasmaRealFloat, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, check, float, PlasmaRealFloat, LDB, N, NRHS );
    PLASMA_splgsy_Tile( (float)N, descA, 51 );

    /* Save A for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, float, LDA, N );

    /* PLASMA SPOSV */
    START_TIMING();
    PLASMA_spotrf_Tile(uplo, descA);
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        PLASMA_splrnt_Tile( descB, 7672 );
        PASTE_TILE_TO_LAPACK( descB, B, check, float, LDB, NRHS );

        PLASMA_spotrs_Tile( uplo, descA, descB );

        PASTE_TILE_TO_LAPACK( descB, X, check, float, LDB, NRHS );

        dparam[IPARAM_RES] = s_check_solution(N, N, NRHS, A, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));

        PASTE_CODE_FREE_MATRIX( descB );
        free( A );
        free( B );
        free( X );
      }

    PASTE_CODE_FREE_MATRIX( descA );

    return 0;
}
