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
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, PLASMA_Complex32_t, LDA, N );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, min(M,N), 1 );
    
    /* Initialize Data */
    PLASMA_cplrnt(M, N, A, LDA, 3456);

    /* Save AT in lapack layout for check */
    PASTE_CODE_ALLOCATE_COPY( Acpy, check, PLASMA_Complex32_t, A, LDA, N );
    
    START_TIMING();
    PLASMA_cgetrf( M, N, A, LDA, piv );
    STOP_TIMING();
    
    /* Check the solution */
    if ( check )
    {
        PASTE_CODE_ALLOCATE_MATRIX( X, 1, PLASMA_Complex32_t, LDB, NRHS );
        PLASMA_cplrnt( N, NRHS, X, LDB, 5673 );
        PASTE_CODE_ALLOCATE_COPY( B, 1, PLASMA_Complex32_t, X, LDB, NRHS );
        
        PLASMA_cgetrs( PlasmaNoTrans, N, NRHS, A, LDA, piv, X, LDB );
        
        dparam[IPARAM_RES] = c_check_solution(M, N, NRHS, Acpy, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        
        free( Acpy ); free( B ); free( X );
    }
    
    free( A );
    free( piv );

    return 0;
}
