/**
 *
 * @generated s Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sgemm_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEMM(M, N, K)
#define _FADDS FADDS_GEMM(M, N, K)

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    float alpha, beta;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    LDB = max(K, iparam[IPARAM_LDB]);
    LDC = max(M, iparam[IPARAM_LDC]);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, float, PlasmaRealFloat, LDA, M, K );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, float, PlasmaRealFloat, LDB, K, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descC, 1, float, PlasmaRealFloat, LDC, M, N );

    /* Initialiaze Data */
    PLASMA_splrnt_Tile( descA, 5373 );
    PLASMA_splrnt_Tile( descB, 7672 );
    PLASMA_splrnt_Tile( descC, 6387 );
    
    LAPACKE_slarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_slarnv_work(1, ISEED, 1, &beta);
    
    /* Save C for check */
    PASTE_TILE_TO_LAPACK( descC, C2, check, float, LDC, N );

    START_TIMING();
    PLASMA_sgemm_Tile( PlasmaNoTrans, PlasmaNoTrans, alpha, descA, descB, beta, descC );
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        PASTE_TILE_TO_LAPACK( descA, A, check, float, LDA, K );
        PASTE_TILE_TO_LAPACK( descB, B, check, float, LDB, N );
        PASTE_TILE_TO_LAPACK( descC, C, check, float, LDC, N );

        dparam[IPARAM_RES] = s_check_gemm( PlasmaNoTrans, PlasmaNoTrans, M, N, K,
                                           alpha, A, LDA, B, LDB, beta, C, C2, LDC,
                                           &(dparam[IPARAM_ANORM]), 
                                           &(dparam[IPARAM_BNORM]), 
                                           &(dparam[IPARAM_XNORM]));
        free(A); free(B); free(C); free(C2);
    }

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );
    PASTE_CODE_FREE_MATRIX( descC );
    return 0;
}
