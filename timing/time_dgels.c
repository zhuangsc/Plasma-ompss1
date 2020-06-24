/**
 *
 * @generated d Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgels"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GEQRF( M, N ) + FMULS_GEQRS( M, N, NRHS ))
#define _FADDS (FADDS_GEQRF( M, N ) + FADDS_GEQRS( M, N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_desc *T;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A,    1,     double, LDA, N   );
    PASTE_CODE_ALLOCATE_MATRIX( x,    1,     double, LDB, NRHS);
    PASTE_CODE_ALLOCATE_MATRIX( Acpy, check, double, LDA, N   );
    PASTE_CODE_ALLOCATE_MATRIX( b,    check, double, LDB, NRHS);

     /* Initialiaze Data */
    PLASMA_dplrnt( M, N,    A, LDA,  453 );
    PLASMA_dplrnt( M, NRHS, x, LDB, 5673 );

    PLASMA_Alloc_Workspace_dgels(M, N, &T);

    /* Save A and b  */
    if (check) {
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'A', M, N,    A, LDA, Acpy, LDA);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'A', M, NRHS, x, LDB, b,    LDB);
    }

    START_TIMING();
    PLASMA_dgels( PlasmaNoTrans, M, N, NRHS, A, LDA, T, x, LDB );
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = d_check_solution(M, N, NRHS, Acpy, LDA, b, x, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(Acpy); free(b);
    }

    PLASMA_Dealloc_Handle_Tile( &T );
    free( A );
    free( x );

    return 0;
}
