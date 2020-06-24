/**
 *
 * @generated c Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cheevd_Tile"
/* See Lawn 41 page 120 */
#define _FMULS ((2. / 3.) * ((float)N * (float)N * (float)N))
#define _FADDS ((2. / 3.) * ((float)N * (float)N * (float)N))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_desc *descT;
    int LDQ;
    int INFO;
    float abstol = LAPACKE_slamch_work('s'); /* not used in this version */
    int nbcomputedeig = 0;
    float vl = 0.0;
    float vu = 0.0;
    int il    = 0;
    int iu    = (int)(0.2*(float)N)/N; // this will compute the first 20% of the eigenvectors. note that you will have to switch range  = PlasmaIvec;
    PLASMA_enum range  = PlasmaAllVec;
    PLASMA_enum uplo   = PlasmaLower;
    PLASMA_enum vec    = PlasmaVec;


    LDA = max(LDA, N);
    LDQ = N;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex32_t, PlasmaComplexFloat, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( Q, (vec == PlasmaVec), PLASMA_Complex32_t, LDA, N );
    PASTE_CODE_ALLOCATE_MATRIX( W, 1, float, N, 1 );

    /* Set Q to 0. */
    if (vec == PlasmaVec) {
        memset(Q,0,LDQ*N*sizeof(PLASMA_Complex32_t));
    }

    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_cplghe_Tile((float)0.0, descA, 51 );

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_cheevr(N, N, &descT);

    START_TIMING();
    INFO = PLASMA_cheevr_Tile(vec, range, uplo, descA, vl, vu, il, iu, abstol, &nbcomputedeig, W, descT, Q, LDQ);
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
