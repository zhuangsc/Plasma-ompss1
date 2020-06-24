/**
 *
 * @generated d Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dsyev_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEBRD( M, N )
#define _FADDS FADDS_GEBRD( M, N )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_desc *descT;
    int jobu  = PlasmaNoVec;
    int jobvt = PlasmaNoVec;
    int INFO;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, double, PlasmaRealDouble, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX( VT, (jobvt == PlasmaVec), double, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( U, (jobu == PlasmaVec), double, M, M );
    PASTE_CODE_ALLOCATE_MATRIX( S, 1, double, N, 1 );

    /* Initialiaze Data */
    PLASMA_dplrnt_Tile(descA, 51 );

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_dgesvd(N, N, &descT);

    if ( jobu == PlasmaVec ) {
        LAPACKE_dlaset_work(LAPACK_COL_MAJOR, 'A', M, M, 0., 1., U, M);
    }
    if ( jobvt == PlasmaVec ) {
        LAPACKE_dlaset_work(LAPACK_COL_MAJOR, 'A', N, N, 0., 1., VT, N);
    }


    START_TIMING(); 
    INFO = PLASMA_dgesvd_Tile(jobu, jobvt, descA, S, descT, U, M, VT, N);
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
