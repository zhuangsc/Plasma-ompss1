/**
 *
 * @generated s Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  float
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_sgetrf_reclap"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(M, NRHS)
#define _FADDS FADDS_GETRF(M, NRHS)

#include "../control/common.h"
#include "./timing.c"

void CORE_sgetrf_reclap_init(void);
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
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, float, LDA, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX( ipiv, 1, int, max(M, NRHS), 1 );

    /* Initialiaze Data */
    PLASMA_splrnt(M, NRHS, A, LDA, 3456);

    /* Save A in lapack layout for check */
    PASTE_CODE_ALLOCATE_COPY(   A2,  check, float, A, LDA, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX( ipiv2, check, int, max(M, NRHS), 1 );
    if ( check ) {
        LAPACKE_sgetrf_work(LAPACK_COL_MAJOR, M, NRHS, A2, LDA, ipiv2 );
    }

    plasma = plasma_context_self();
    PLASMA_Sequence_Create(&sequence);
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    QUARK_Task_Flag_Set(&task_flags, TASK_THREAD_COUNT, iparam[IPARAM_THRDNBR] );

    plasma_dynamic_spawn();
    CORE_sgetrf_reclap_init();

    START_TIMING();
    QUARK_CORE_sgetrf_reclap(plasma->quark, &task_flags,
                             M, NRHS, NRHS,
                             A, LDA, ipiv,
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

        /* Check ipiv */
        for(i=0; i<NRHS; i++)
        {
            if( ipiv[i] != ipiv2[i] ) {
                fprintf(stderr, "\nPLASMA (ipiv[%ld] = %d, A[%ld] = %e) / LAPACK (ipiv[%ld] = %d, A[%ld] = [%e])\n",
                        i, ipiv[i],  i, (A[  i * LDA + i ]), 
                        i, ipiv2[i], i, (A2[ i * LDA + i ])); 
                break;
            }
        }

        dparam[IPARAM_ANORM] = LAPACKE_slange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                   M, NRHS, A, LDA, work);
        dparam[IPARAM_XNORM] = LAPACKE_slange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                   M, NRHS, A2, LDA, work);
        dparam[IPARAM_BNORM] = 0.0;

        CORE_sgeadd( M, NRHS, -1.0, A, LDA, A2, LDA);

        dparam[IPARAM_RES] = LAPACKE_slange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaMaxNorm), 
                                                 M, NRHS, A2, LDA, work);

        free( A2 );
        free( ipiv2 );
        free( work );
    }
    
    free( A  );
    free( ipiv );

    return 0;
}
