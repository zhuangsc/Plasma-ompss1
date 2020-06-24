/**
 *
 * @generated s Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sgetrf_tntpiv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(M, N)
#define _FADDS FADDS_GETRF(M, N)

#include "./timing.c"

#define A(m,n)    BLKADDR(A, float, m, n)

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, float, PlasmaRealFloat, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, min(M, N), 1 );

    PLASMA_splrnt_Tile(descA, 3456);

    /* Save AT in lapack layout for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, float, LDA, N );

    START_TIMING();
    PLASMA_sgetrf_tntpiv_Tile( descA, piv );
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, float, PlasmaRealFloat, LDB, N, NRHS );
        PLASMA_splrnt_Tile( descB, 7732 );
        PASTE_TILE_TO_LAPACK( descB, b, check, float, LDB, NRHS );

        PLASMA_sgetrs_Tile( PlasmaNoTrans, descA, piv, descB );

        PASTE_TILE_TO_LAPACK( descB, x, check, float, LDB, NRHS );
        dparam[IPARAM_RES] = s_check_solution(M, N, NRHS, A, LDA, b, x, LDB,
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
