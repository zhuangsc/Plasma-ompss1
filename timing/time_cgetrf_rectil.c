/**
 *
 * @generated c Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cgetrf_rectil"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(M, NRHS)
#define _FADDS FADDS_GETRF(M, NRHS)

#include "../control/common.h"
#include "./timing.c"

void CORE_cgetrf_rectil_init(void);
extern plasma_context_t*  plasma_context_self(void);

/*
 * WARNING: the check is only working with LAPACK Netlib
 * which choose the same pivot than this code.
 * MKL has a different code and can pick a different pivot 
 * if two elments have the same absolute value but not the 
 * same sign for example.
 */

static int
RunTest(int *iparam, float *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    plasma_context_t *plasma;
    Quark_Task_Flags  task_flags = Quark_Task_Flags_Initializer;
    PLASMA_sequence  *sequence = NULL;
    PLASMA_request    request = PLASMA_REQUEST_INITIALIZER;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex32_t, PlasmaComplexFloat, LDA, M, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX( ipiv, 1, int, max(M, NRHS), 1 );

    /* Initialiaze Data */
    PLASMA_cplrnt_Tile(descA, 3456);

    /* Save A in lapack layout for check */
    PASTE_TILE_TO_LAPACK( descA, A2, check, PLASMA_Complex32_t, LDA, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX( ipiv2, check, int, max(M, NRHS), 1 );

    /* Save AT in lapack layout for check */
    if ( check ) {
        LAPACKE_cgetrf_work(LAPACK_COL_MAJOR, M, NRHS, A2, LDA, ipiv2 );
    }

    plasma = plasma_context_self();
    PLASMA_Sequence_Create(&sequence);
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    QUARK_Task_Flag_Set(&task_flags, TASK_THREAD_COUNT, iparam[IPARAM_THRDNBR] );

    plasma_dynamic_spawn();
    CORE_cgetrf_rectil_init();

    START_TIMING();
    QUARK_CORE_cgetrf_rectil(plasma->quark, &task_flags,
                             *descA, descA->mat, descA->mb*descA->nb, ipiv,
                             sequence, &request,
                             0, 0,
                             iparam[IPARAM_THRDNBR]);
    PLASMA_Sequence_Wait(sequence);
    STOP_TIMING();
    
    PLASMA_Sequence_Destroy(sequence);

    /* Check the solution */
    if ( check )
    {
        int64_t i;
        float *work = (float *)malloc(max(M, NRHS)*sizeof(float));
        PASTE_TILE_TO_LAPACK( descA, A, 1, PLASMA_Complex32_t, LDA, NRHS );

        /* Check ipiv */
        for(i=0; i<NRHS; i++)
        {
            if( ipiv[i] != ipiv2[i] ) {
                fprintf(stderr, "\nPLASMA (ipiv[%ld] = %d, A[%ld] = %e) / LAPACK (ipiv[%ld] = %d, A[%ld] = [%e])\n",
                        i, ipiv[i],  i, crealf(A[  i * LDA + i ]), 
                        i, ipiv2[i], i, crealf(A2[ i * LDA + i ])); 
                break;
            }
        }

        dparam[IPARAM_ANORM] = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                   M, NRHS, A, LDA, work);
        dparam[IPARAM_XNORM] = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                   M, NRHS, A2, LDA, work);
        dparam[IPARAM_BNORM] = 0.0;

        CORE_cgeadd( M, NRHS, -1.0, A, LDA, A2, LDA);

        dparam[IPARAM_RES] = LAPACKE_clange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                 M, NRHS, A2, LDA, work);

        free( A );
        free( A2 );
        free( ipiv2 );
        free( work );
    }
    
    /* Deallocate Workspace */
    PASTE_CODE_FREE_MATRIX( descA );
    free( ipiv );

    return 0;
}
