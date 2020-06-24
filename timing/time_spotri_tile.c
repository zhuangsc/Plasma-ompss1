/**
 *
 * @generated s Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_spotri_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( N ) + FMULS_POTRI( N ))
#define _FADDS (FADDS_POTRF( N ) + FADDS_POTRI( N ))

//#define POTRI_SYNC

#include "./timing.c"

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_enum uplo = PlasmaLower;

    LDA = max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, float, PlasmaRealFloat, LDA, N, N );

    /* 
     * Initialize Data 
     * It's done in static to avoid having the same sequence than one 
     * the function we want to trace
     */
    PLASMA_splgsy_Tile( (float)N, descA, 51 );
    
    /* PLASMA SPOTRF / STRTRI / SLAUUM  */
    /*
     * Example of the different way to combine several asynchonous calls
     */
#if defined(TRACE_BY_SEQUENCE)
    {
        PLASMA_sequence *sequence[3];
        PLASMA_request request[3] = { PLASMA_REQUEST_INITIALIZER, 
                                      PLASMA_REQUEST_INITIALIZER, 
                                      PLASMA_REQUEST_INITIALIZER };
        
        PLASMA_Sequence_Create(&sequence[0]);
        PLASMA_Sequence_Create(&sequence[1]);
        PLASMA_Sequence_Create(&sequence[2]);
        
        if ( ! iparam[IPARAM_ASYNC] ) {
            START_TIMING();

            PLASMA_spotrf_Tile_Async(uplo, descA,                sequence[0], &request[0]);
            PLASMA_Sequence_Wait(sequence[0]);
            
            PLASMA_strtri_Tile_Async(uplo, PlasmaNonUnit, descA, sequence[1], &request[1]);
            PLASMA_Sequence_Wait(sequence[1]);
            
            PLASMA_slauum_Tile_Async(uplo, descA,                sequence[2], &request[2]);
            PLASMA_Sequence_Wait(sequence[2]);
            STOP_TIMING();

        } else {

            START_TIMING();
            PLASMA_spotrf_Tile_Async(uplo, descA,                sequence[0], &request[0]);
            PLASMA_strtri_Tile_Async(uplo, PlasmaNonUnit, descA, sequence[1], &request[1]);
            PLASMA_slauum_Tile_Async(uplo, descA,                sequence[2], &request[2]);

            PLASMA_Sequence_Wait(sequence[0]);
            PLASMA_Sequence_Wait(sequence[1]);
            PLASMA_Sequence_Wait(sequence[2]);
            STOP_TIMING();
        }
        
        PLASMA_Sequence_Destroy(sequence[0]);
        PLASMA_Sequence_Destroy(sequence[1]);
        PLASMA_Sequence_Destroy(sequence[2]);
    }       
#else
    {
        if ( ! iparam[IPARAM_ASYNC] ) {

            START_TIMING();
            PLASMA_spotrf_Tile(uplo, descA);
            PLASMA_strtri_Tile(uplo, PlasmaNonUnit, descA);
            PLASMA_slauum_Tile(uplo, descA);
            STOP_TIMING();

        } else {

            /* Default: we use Asynchonous call with only one sequence */
            PLASMA_sequence *sequence;
            PLASMA_request request[2] = { PLASMA_REQUEST_INITIALIZER, 
                                          PLASMA_REQUEST_INITIALIZER };
        
            START_TIMING();
            PLASMA_Sequence_Create(&sequence);
            PLASMA_spotrf_Tile_Async(uplo, descA, sequence, &request[0]);
            PLASMA_spotri_Tile_Async(uplo, descA, sequence, &request[1]);
            PLASMA_Sequence_Wait(sequence);
            STOP_TIMING();
        
            PLASMA_Sequence_Destroy(sequence);       
        }
    }
#endif
    
    /* Check the solution */
    if ( check )
    {
        dparam[IPARAM_ANORM] = 0.0;
        dparam[IPARAM_XNORM] = 0.0;
        dparam[IPARAM_BNORM] = 0.0;
        dparam[IPARAM_RES]   = 0.0;
    }

    PASTE_CODE_FREE_MATRIX( descA );

    return 0;
}
