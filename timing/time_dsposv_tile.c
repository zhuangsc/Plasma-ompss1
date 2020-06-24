/**
 *
 * @generated ds Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dposv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( N ) + FMULS_POTRS( N, NRHS ))
#define _FADDS (FADDS_POTRF( N ) + FADDS_POTRS( N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    int iter;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, double, PlasmaRealDouble, LDA, N, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, double, PlasmaRealDouble, LDB, N, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descX, 1, double, PlasmaRealDouble, LDB, N, NRHS );

    PLASMA_dplgsy_Tile((double)N, descA, 51 );
    PLASMA_dplrnt_Tile(descB, 5723);

    /* Save AT and bT in lapack layout for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, double, LDA, N    );
    PASTE_TILE_TO_LAPACK( descB, b, check, double, LDB, NRHS );

    /* PLASMA DSPOSV */
    START_TIMING();
    PLASMA_dsposv_Tile(PlasmaUpper, descA, descB, descX, &iter);
    STOP_TIMING();

    /* Check the solution */
    if (check)
    {
        PASTE_TILE_TO_LAPACK( descX, x, check, double, LDB, NRHS );

        dparam[IPARAM_RES] = d_check_solution(N, N, NRHS, A, LDA, b, x, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(A); free(b); free(x);
    }

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );
    PASTE_CODE_FREE_MATRIX( descX );

    return 0;
}
