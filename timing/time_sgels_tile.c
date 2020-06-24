/**
 *
 * @generated s Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sgels_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GEQRF( M, N ) + FMULS_GEQRS( M, N, NRHS ))
#define _FADDS (FADDS_GEQRF( M, N ) + FADDS_GEQRS( M, N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PLASMA_desc *descT;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, float, PlasmaRealFloat, LDA, M, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, float, PlasmaRealFloat, LDB, M, NRHS );

    PLASMA_splrnt_Tile( descA, 5373 );
    PLASMA_splrnt_Tile( descB,  673 );

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_sgels_Tile(M, N, &descT);

    /* Save A and B for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, float, LDA, N    );
    PASTE_TILE_TO_LAPACK( descB, B, check, float, LDB, NRHS );

    /* Do the computations */
    START_TIMING();
    PLASMA_sgels_Tile( PlasmaNoTrans, descA, descT, descB );
    STOP_TIMING();
    
    /* Allocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    /* Check the solution */
    if ( check )
      {
        /* Copy solution to X */
        PASTE_TILE_TO_LAPACK( descB, X, 1, float, LDB, NRHS );

        dparam[IPARAM_RES] = s_check_solution(M, N, NRHS, A, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(A); free(B); free(X);
      }

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );

    return 0;
}
