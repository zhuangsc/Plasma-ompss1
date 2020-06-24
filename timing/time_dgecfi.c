/**
 *
 * @generated d Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgecfi"
/* See Lawn 41 page 120 */
#define _FMULS (0.0)
#define _FADDS (M * N * sizeof(double))

#include "./timing.c"

int d_check_conversion(int m, int n, int mba, int nba, int mbb, int nbb,
                      double *A, double *B, 
                      int (*mapA)(int, int, int, int, int, int), int (*mapB)(int, int, int, int, int, int)) {
    int i, j;

    for( j=0; j<n; j++) {
        for (i=0; i<m; i++) {
            if (A[ mapA(m, n, mba, nba, i, j) ] != B[ mapB(m, n, mbb, nbb, i, j) ] ) {
                return -1; 
            }
        }
    }
    return 0;
}

static int
RunTest(int *iparam, _PREC *dparam, real_Double_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    LDA = M;

    dparam[IPARAM_ANORM] = (_PREC)M;
    dparam[IPARAM_BNORM] = (_PREC)_FADDS;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, double, LDA, N );
    
    /* Initialize Data */
    PLASMA_dplrnt(M, N, A, LDA, 3456);

    /* Save AT in lapack layout for check */
    PASTE_CODE_ALLOCATE_COPY( Acpy, check, double, A, LDA, N );
    
    START_TIMING();
    PLASMA_dgecfi( M, N, A, PlasmaCM, M, 1, PlasmaCCRB, MB, NB);
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = (_PREC)d_check_conversion(M, N, M, 1, MB, NB, Acpy, A, map_CM, map_CCRB);
        free(Acpy);
    }

    free( A );
    return 0;
}
