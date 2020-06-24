/**
 *
 * @generated c Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_chegv_Tile"
/* See Lawn 41 page 120 */
/* cholesky + 2 trsm's + trd */
#define _FMULS (((2. / 3.) * ((float)N * (float)N * (float)N)) + (N * (1.0 / 6.0 * N + 0.5) * N) + 2 * (N * N * (float)((N + 1) / 2.0) ))
#define _FADDS (((2. / 3.) * ((float)N * (float)N * (float)N)) + (N * (1.0 / 6.0 * N )      * N) + 2 * (N * N * (float)((N + 1) / 2.0) ))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_desc *descT;
    int itype = 1;
    int vec   = PlasmaNoVec;
    int uplo  = PlasmaUpper;

    LDA = max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex32_t, PlasmaComplexFloat, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, PLASMA_Complex32_t, PlasmaComplexFloat, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descQ, (vec == PlasmaVec), PLASMA_Complex32_t, PlasmaComplexFloat, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( W, 1, float, N, 1 );

    /* Initialization */
    PLASMA_cplghe_Tile((float)0.0, descA, 51   );
    PLASMA_cplghe_Tile((float)N,   descB, 3753 );

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_chegv(N, N, &descT);

    START_TIMING();
    PLASMA_chegv_Tile( itype, vec, uplo, descA, descB, W, descT, descQ );
    STOP_TIMING();

    /* DeAllocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    if (vec == PlasmaVec) {
        PASTE_CODE_FREE_MATRIX( descQ );
    }
    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );
    free( W );

    return 0;
}
