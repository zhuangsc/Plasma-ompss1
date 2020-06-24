/**
 *
 * @generated c Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cposv_Tile"
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
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex32_t, PlasmaComplexFloat, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, PLASMA_Complex32_t, PlasmaComplexFloat, LDB, N, NRHS );

    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_cplghe_Tile((float)N, descA, 51 );
    PLASMA_cplrnt_Tile( descB, 7732 );

    /* Save AT and bT in lapack layout for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, PLASMA_Complex32_t, LDA, N );
    PASTE_TILE_TO_LAPACK( descB, B, check, PLASMA_Complex32_t, LDB, NRHS );

    /* PLASMA CPOSV */
    START_TIMING();
    PLASMA_cposv_Tile(uplo, descA, descB);
    STOP_TIMING();

    /* Check the solution */
    if (check)
      {
        PASTE_TILE_TO_LAPACK( descB, X, check, PLASMA_Complex32_t, LDB, NRHS );

        dparam[IPARAM_RES] = c_check_solution(N, N, NRHS, A, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(A); free(B); free(X);
      }

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );

    return 0;
}
