/**
 *
 * @generated c Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cgetrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(M, N)
#define _FADDS FADDS_GETRF(M, N)

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex32_t, PlasmaComplexFloat, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, min(M, N), 1 );

    PLASMA_cplrnt_Tile(descA, 3456);

    /* Save AT in lapack layout for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, PLASMA_Complex32_t, LDA, N );

    START_TIMING();
    PLASMA_cgetrf_Tile( descA, piv );
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, PLASMA_Complex32_t, PlasmaComplexFloat, LDB, N, NRHS );
        PLASMA_cplrnt_Tile( descB, 7732 );
        PASTE_TILE_TO_LAPACK( descB, b, check, PLASMA_Complex32_t, LDB, NRHS );

        PLASMA_cgetrs_Tile( PlasmaNoTrans, descA, piv, descB );

        PASTE_TILE_TO_LAPACK( descB, x, check, PLASMA_Complex32_t, LDB, NRHS );
        dparam[IPARAM_RES] = c_check_solution(M, N, NRHS, A, LDA, b, x, LDB,
                                              &(dparam[IPARAM_ANORM]),
                                              &(dparam[IPARAM_BNORM]),
                                              &(dparam[IPARAM_XNORM]));

        PASTE_CODE_FREE_MATRIX( descB );
        free(A); free(b); free(x);
    }

    PASTE_CODE_FREE_MATRIX( descA );
    free( piv );

    return 0;
}
