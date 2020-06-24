/**
 *
 * @generated c Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cgeqrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEQRF( M, N )
#define _FADDS FADDS_GEQRF( M, N )

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PLASMA_desc *descT;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex32_t, PlasmaComplexFloat, LDA, M, N );
    PLASMA_cplrnt_Tile( descA, 5373 );

    /* Save A for check */
    PASTE_TILE_TO_LAPACK( descA, A, ( check && M == N ), PLASMA_Complex32_t, LDA, N );
    
    /* Allocate B for check */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, (check && M == N), PLASMA_Complex32_t, PlasmaComplexFloat, LDB, M, NRHS );
     
    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_cgels_Tile(M, N, &descT);

    /* Do the computations */
    START_TIMING();
    PLASMA_cgeqrf_Tile( descA, descT );
    STOP_TIMING();
    
    /* Check the solution */
    if ( check && M == N )
    {
        /* Initialize and save B */
        PLASMA_cplrnt_Tile( descB, 2264 );
        PASTE_TILE_TO_LAPACK( descB, B, 1, PLASMA_Complex32_t, LDB, NRHS );

        /* Compute the solution */
        PLASMA_cgeqrs_Tile( descA, descT, descB );

        /* Copy solution to X */
        PASTE_TILE_TO_LAPACK( descB, X, 1, PLASMA_Complex32_t, LDB, NRHS );

        dparam[IPARAM_RES] = c_check_solution(M, N, NRHS, A, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));

        /* Free checking structures */
        PASTE_CODE_FREE_MATRIX( descB );
        free( A ); 
        free( B ); 
        free( X );
    }

    /* Free data */
    PLASMA_Dealloc_Handle_Tile(&descT);
    PASTE_CODE_FREE_MATRIX( descA );

    return 0;
}
