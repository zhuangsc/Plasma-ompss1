
/**
 *
 * @generated d Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dsyevd_Tile"
/* See Lawn 41 page 120 */
#define _FMULS ((2. / 3.) * ((double)N * (double)N * (double)N))
#define _FADDS ((2. / 3.) * ((double)N * (double)N * (double)N))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_desc *descT;
   // double *Q;
    int LDQ;
    int uplo = PlasmaLower;
    int vec  = PlasmaVec;
    int INFO;

    LDA = max(LDA, N);
    LDQ = N;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, double, PlasmaRealDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( Q, (vec == PlasmaVec), double, LDA, N );
    PASTE_CODE_ALLOCATE_MATRIX( W, 1, double, N, 1 );

    /* Set Q to 0. */
    if (vec == PlasmaVec) {
        memset(Q,0,LDQ*N*sizeof(double));
    }

    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_dplgsy_Tile((double)0.0, descA, 51 );
       
    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_dsyevd(N, N, &descT);

    START_TIMING();
    INFO = PLASMA_dsyevd_Tile( vec, uplo, descA, W, descT, Q, LDQ );
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
    free( W );

    return 0;
}
