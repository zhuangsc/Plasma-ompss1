/**
 *
 * @generated d Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgesv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRF( N, N ) + FMULS_GETRS( N, NRHS ))
#define _FADDS (FADDS_GETRF( N, N ) + FADDS_GETRS( N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }
    
    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, double, PlasmaRealDouble, LDA, N, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, double, PlasmaRealDouble, LDB, N, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, N, 1 );
    
    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_dplrnt_Tile( descA, 8796 );
    PLASMA_dplrnt_Tile( descB, 7732 );

    /* Save AT and bT in lapack layout for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, double, LDA, N    );
    PASTE_TILE_TO_LAPACK( descB, b, check, double, LDB, NRHS );

    START_TIMING();
    PLASMA_dgesv_Tile( descA, piv, descB );
    STOP_TIMING();
    
    /* Check the solution */
    if ( check )
    {
        PASTE_TILE_TO_LAPACK( descB, x, check, double, LDB, NRHS );
        
        dparam[IPARAM_RES] = d_check_solution(N, N, NRHS, A, LDA, b, x, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(A); free(b); free(x);
    }

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );
    free( piv );

    return 0;
}
