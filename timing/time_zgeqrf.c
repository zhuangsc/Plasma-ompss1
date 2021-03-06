/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgeqrf"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEQRF(M, N)
#define _FADDS FADDS_GEQRF(M, N)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_desc *T;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, PLASMA_Complex64_t, LDA, N );

    /* Initialize Data */
    PLASMA_zplrnt(M, N, A, LDA, 3456);

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_zgels(M, N, &T);

    /* Save AT in lapack layout for check */
    PASTE_CODE_ALLOCATE_COPY( Acpy, check, PLASMA_Complex64_t, A, LDA, N );

    START_TIMING();
    PLASMA_zgeqrf( M, N, A, LDA, T );
    STOP_TIMING();
    
    /* Check the solution */
    if ( check )
    {
        PASTE_CODE_ALLOCATE_MATRIX( X, 1, PLASMA_Complex64_t, LDB, NRHS );
        PLASMA_zplrnt( N, NRHS, X, LDB, 5673 );
        PASTE_CODE_ALLOCATE_COPY( B, 1, PLASMA_Complex64_t, X, LDB, NRHS );
        
        PLASMA_zgeqrs(M, N, NRHS, A, LDA, T, X, LDB);

        dparam[IPARAM_RES] = z_check_solution(M, N, NRHS, Acpy, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));

        free( Acpy ); 
        free( B ); 
        free( X );
      }

    /* Free Workspace */
    PLASMA_Dealloc_Handle_Tile( &T );
    free( A );

    return 0;
}
