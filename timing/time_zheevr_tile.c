/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zheevd_Tile"
/* See Lawn 41 page 120 */
#define _FMULS ((2. / 3.) * ((double)N * (double)N * (double)N))
#define _FADDS ((2. / 3.) * ((double)N * (double)N * (double)N))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_desc *descT;
    int LDQ;
    int INFO;
    double abstol = LAPACKE_dlamch_work('s'); /* not used in this version */
    int nbcomputedeig = 0;
    double vl = 0.0;
    double vu = 0.0;
    int il    = 0;
    int iu    = (int)(0.2*(double)N)/N; // this will compute the first 20% of the eigenvectors. note that you will have to switch range  = PlasmaIvec;
    PLASMA_enum range  = PlasmaAllVec;
    PLASMA_enum uplo   = PlasmaLower;
    PLASMA_enum vec    = PlasmaVec;


    LDA = max(LDA, N);
    LDQ = N;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( Q, (vec == PlasmaVec), PLASMA_Complex64_t, LDA, N );
    PASTE_CODE_ALLOCATE_MATRIX( W, 1, double, N, 1 );

    /* Set Q to 0. */
    if (vec == PlasmaVec) {
        memset(Q,0,LDQ*N*sizeof(PLASMA_Complex64_t));
    }

    /* Initialize AT and bT for Symmetric Positif Matrix */
    PLASMA_zplghe_Tile((double)0.0, descA, 51 );

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_zheevr(N, N, &descT);

    START_TIMING();
    INFO = PLASMA_zheevr_Tile(vec, range, uplo, descA, vl, vu, il, iu, abstol, &nbcomputedeig, W, descT, Q, LDQ);
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
