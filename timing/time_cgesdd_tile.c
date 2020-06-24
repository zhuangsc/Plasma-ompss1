/**
 *
 * @generated c Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cheev_Tile"
/* See Lawn 41 page 120 */
#define _FMULS ((2. / 3.) * ((float)N * (float)N * (float)N))
#define _FADDS ((2. / 3.) * ((float)N * (float)N * (float)N))

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_desc *descT;
    int jobu  = PlasmaVec;
    int jobvt = PlasmaVec;
    int INFO;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex32_t, PlasmaComplexFloat, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX( U, (jobu == PlasmaVec), PLASMA_Complex32_t, M, M );	    
    PASTE_CODE_ALLOCATE_MATRIX( VT, (jobvt == PlasmaVec), PLASMA_Complex32_t, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( S, 1, float, N, 1 );

    /* Initialiaze Data */
    PLASMA_cplrnt_Tile(descA, 51 );

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_cgesvd(N, N, &descT);



//printf("\nstart test\n");

    if ( jobu == PlasmaVec ) {
        LAPACKE_claset_work(LAPACK_COL_MAJOR, 'A', M, M, 0., 1., U, M);
    }
    if ( jobvt == PlasmaVec ) {
        LAPACKE_claset_work(LAPACK_COL_MAJOR, 'A', N, N, 0., 1., VT, N);
    }


    START_TIMING(); 
    INFO = PLASMA_cgesdd_Tile(jobu, jobvt, descA, S, descT, U, M, VT, N);
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

    if (jobu == PlasmaVec) {
      free( U );
    }
    if (jobvt == PlasmaVec) {
      free( VT );
    }
    PASTE_CODE_FREE_MATRIX( descA );
    free( S );

    return 0;
}
