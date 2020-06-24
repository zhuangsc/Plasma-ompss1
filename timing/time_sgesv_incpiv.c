/**
 *
 * @generated s Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sgesv_incpiv"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRF( N, N ) + FMULS_GETRS( N, NRHS ))
#define _FADDS (FADDS_GETRF( N, N ) + FADDS_GETRS( N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PLASMA_desc *L;
    int *piv;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }
    
    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, float, LDA, N    );
    PASTE_CODE_ALLOCATE_MATRIX( X, 1, float, LDB, NRHS );
    
    /* Initialiaze Data */
    PLASMA_splrnt( N, N,    A, LDA,   51 );
    PLASMA_splrnt( N, NRHS, X, LDB, 5673 );

    PLASMA_Alloc_Workspace_sgesv_incpiv(N, &L, &piv);

    /* Save A and b  */
    PASTE_CODE_ALLOCATE_COPY( Acpy, check, float, A, LDA, N    );
    PASTE_CODE_ALLOCATE_COPY( B,    check, float, X, LDB, NRHS );

    START_TIMING();
    PLASMA_sgesv_incpiv( N, NRHS, A, LDA, L, piv, X, LDB );
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = s_check_solution(N, N, NRHS, Acpy, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(Acpy); free(B);
    }

    PLASMA_Dealloc_Handle_Tile( &L );
    free( piv );
    free( X );
    free( A );


    return 0;
}
