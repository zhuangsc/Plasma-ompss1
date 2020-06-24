/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgetrf_tntpiv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(M, N)
#define _FADDS FADDS_GETRF(M, N)

#include "./timing.c"

#define A(m,n)    BLKADDR(A, PLASMA_Complex64_t, m, n)

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, min(M, N), 1 );

    PLASMA_zplrnt_Tile(descA, 3456);

    /* Save AT in lapack layout for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, PLASMA_Complex64_t, LDA, N );

    START_TIMING();
    PLASMA_zgetrf_tntpiv_Tile( descA, piv );
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDB, N, NRHS );
        PLASMA_zplrnt_Tile( descB, 7732 );
        PASTE_TILE_TO_LAPACK( descB, b, check, PLASMA_Complex64_t, LDB, NRHS );

        PLASMA_zgetrs_Tile( PlasmaNoTrans, descA, piv, descB );

        PASTE_TILE_TO_LAPACK( descB, x, check, PLASMA_Complex64_t, LDB, NRHS );
        dparam[IPARAM_RES] = z_check_solution(M, N, NRHS, A, LDA, b, x, LDB,
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
