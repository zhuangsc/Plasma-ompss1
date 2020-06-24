/**
 *
 * @generated c Tue Jan  7 11:45:23 2014
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_ctrsm"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_TRSM( PlasmaLeft, N, NRHS )
#define _FADDS FADDS_TRSM( PlasmaLeft, N, NRHS )

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_)
{
    PLASMA_Complex32_t alpha;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    LDA = max( LDA, N );

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A,      1, PLASMA_Complex32_t, LDA, N   );
    PASTE_CODE_ALLOCATE_MATRIX( B,      1, PLASMA_Complex32_t, LDB, NRHS);
    PASTE_CODE_ALLOCATE_MATRIX( B2, check, PLASMA_Complex32_t, LDB, NRHS);

     /* Initialiaze Data */
    PLASMA_cplgsy( (PLASMA_Complex32_t)N, N, A, LDA, 453 );
    PLASMA_cplrnt( N, NRHS, B, LDB, 5673 );
    LAPACKE_clarnv_work(1, ISEED, 1, &alpha);
    alpha = 10.; /*alpha * N  /  2.;*/

    if (check)
    {
        memcpy(B2, B, LDB*NRHS*sizeof(PLASMA_Complex32_t));
    }

    START_TIMING();
    PLASMA_ctrsm( PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaUnit,
                  N, NRHS, alpha, A, LDA, B, LDB );
    STOP_TIMING();

    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = c_check_trsm( PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaUnit, 
                                           N, NRHS,
                                           alpha, A, LDA, B, B2, LDB,
                                           &(dparam[IPARAM_ANORM]), 
                                           &(dparam[IPARAM_BNORM]),
                                           &(dparam[IPARAM_XNORM]));
        free(B2);
    }

    free( A );
    free( B );

    return 0;
}
