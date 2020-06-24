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

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, double, LDA, N );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, min(M,N), 1 );
    
    /* Initialize Data */
    PLASMA_dplrnt(M, N, A, LDA, 3456);

    /* Save AT in lapack layout for check */
    PASTE_CODE_ALLOCATE_COPY( Acpy, check, double, A, LDA, N );
    
    START_TIMING();
    PLASMA_dgetrf( M, N, A, LDA, piv );
    STOP_TIMING();
    
    /* Check the solution */
    if ( check )
    {
        PASTE_CODE_ALLOCATE_MATRIX( X, 1, double, LDB, NRHS );
        PLASMA_dplrnt( N, NRHS, X, LDB, 5673 );
        PASTE_CODE_ALLOCATE_COPY( B, 1, double, X, LDB, NRHS );
        
        PLASMA_dgetrs( PlasmaNoTrans, N, NRHS, A, LDA, piv, X, LDB );
        
        dparam[IPARAM_RES] = d_check_solution(M, N, NRHS, Acpy, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        
        free( Acpy ); free( B ); free( X );
    }
    
    free( A );
    free( piv );

    return 0;
}
