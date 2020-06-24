
/**
 *
 * @generated s Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_ssytrd_Tile"
/* See Lawn 41 page 120 */
#define _FMULS ((2. / 3.) * ((float)N * (float)N * (float)N))
#define _FADDS ((2. / 3.) * ((float)N * (float)N * (float)N))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_desc *descT;
   // float *Q;
    int LDQ;
    int uplo = PlasmaLower;
    int vec  = PlasmaVec;
    int INFO;

    LDA = max(LDA, N);
    LDQ = N;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, float, PlasmaRealFloat, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( Q, (vec == PlasmaVec), float, LDA, N );
    PASTE_CODE_ALLOCATE_MATRIX( D, 1, float, N, 1 );
    PASTE_CODE_ALLOCATE_MATRIX( E, 1, float, N, 1 );

    /* Set Q to 0. */
    if (vec == PlasmaVec) {
        memset(Q,0,LDQ*N*sizeof(float));
    }

    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_splgsy_Tile((float)0.0, descA, 51 );
       
    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_ssyevd(N, N, &descT);

    START_TIMING();
    INFO = PLASMA_ssytrd_Tile( vec, uplo, descA, D, E, descT, Q, LDQ );
    STOP_TIMING();
    if(INFO!=0){
            printf(" ERROR OCCURED INFO %d\n",INFO);
    }
 
    /* Check the solution */
    if ( check )
      {
      }

    /* DeAllocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    if (vec == PlasmaVec) {
      free( Q );
    }
    PASTE_CODE_FREE_MATRIX( descA );
    free( D );
    free( E );

    return 0;
}
