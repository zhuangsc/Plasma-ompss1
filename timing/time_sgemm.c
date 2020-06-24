/**
 *
 * @generated s Tue Jan  7 11:45:23 2014
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sgemm"
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
    PASTE_CODE_ALLOCATE_MATRIX( A,      1, float, LDA, K   );
    PASTE_CODE_ALLOCATE_MATRIX( B,      1, float, LDB, N   );
    PASTE_CODE_ALLOCATE_MATRIX( C,      1, float, LDC, N   );
    PASTE_CODE_ALLOCATE_MATRIX( C2, check, float, LDC, N   );

    PLASMA_splrnt( M, K, A, LDA,  453 );
    PLASMA_splrnt( K, N, B, LDB, 5673 );
    PLASMA_splrnt( M, N, C, LDC,  740 );

    LAPACKE_slarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_slarnv_work(1, ISEED, 1, &beta );

    if (check)
    {
        memcpy(C2, C, LDC*N*sizeof(float));
    }

    START_TIMING();
    PLASMA_sgemm( PlasmaNoTrans, PlasmaNoTrans, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC );
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = s_check_gemm( PlasmaNoTrans, PlasmaNoTrans, M, N, K,
                                           alpha, A, LDA, B, LDB, beta, C, C2, LDC,
                                           &(dparam[IPARAM_ANORM]), 
                                           &(dparam[IPARAM_BNORM]), 
                                           &(dparam[IPARAM_XNORM]));
        free(C2);
    }

    free( A );
    free( B );
    free( C );

    return 0;
}
