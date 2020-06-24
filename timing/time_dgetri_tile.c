/**
 *
 * @generated d Tue Jan  7 11:45:24 2014
 *
 **/
#define _TYPE  double
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_dgetri_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRF(M, N) + FMULS_GETRI( N ))
#define _FADDS (FADDS_GETRF(M, N) + FADDS_GETRI( N ))

//#define GETRI_SYNC

#include "./timing.c"

/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */
#if 0
static int check_getri_factorization(PLASMA_desc *descA1, PLASMA_desc *descA2, int *IPIV)
{
    int info_factorization;
    double Rnorm, Anorm, Xnorm, Bnorm, result;
    double eps = LAPACKE_dlamch_work('e');
    PLASMA_desc        *descB, *descX;
    double *b = (double *)malloc((descA1->m)*sizeof(double));
    double *x = (double *)malloc((descA1->m)*sizeof(double));

    PLASMA_Desc_Create(&descB, b, PlasmaRealDouble, descA1->mb, descA1->nb, descA1->bsiz,
                       descA1->m, 1, 0, 0, descA1->m, 1);
    PLASMA_Desc_Create(&descX, x, PlasmaRealDouble, descA1->mb, descA1->nb, descA1->bsiz,
                       descA1->m, 1, 0, 0, descA1->m, 1);

    PLASMA_dplrnt_Tile( descX, 537 );
    PLASMA_dlacpy_Tile( PlasmaUpperLower, descX, descB);

    PLASMA_dgetrs_Tile( PlasmaNoTrans, descA2, IPIV, descX );

    Xnorm = PLASMA_dlange_Tile(PlasmaInfNorm, descX);
    Anorm = PLASMA_dlange_Tile(PlasmaInfNorm, descA1);
    Bnorm = PLASMA_dlange_Tile(PlasmaInfNorm, descB);

    PLASMA_dgemm_Tile( PlasmaNoTrans, PlasmaNoTrans,
                       (double)1.,  descA1, descX,
                       (double)-1., descB);

    Rnorm = PLASMA_dlange_Tile(PlasmaInfNorm, descB);

    if (getenv("PLASMA_TESTING_VERBOSE"))
      printf( "||A||_oo=%f\n||X||_oo=%f\n||B||_oo=%f\n||A X - B||_oo=%e\n", Anorm, Xnorm, Bnorm, Rnorm );

    result = Rnorm / ( (Anorm*Xnorm+Bnorm)*(descA1->m)*eps ) ;
    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);

    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- The factorization is suspicious ! \n");
        info_factorization = 1;
     }
    else{
        printf("-- The factorization is CORRECT ! \n");
        info_factorization = 0;
    }
    free(x); free(b);
    PLASMA_Desc_Destroy(&descB);
    PLASMA_Desc_Destroy(&descX);

    return info_factorization;
}
#endif

/*------------------------------------------------------------------------
 *  Check the accuracy of the computed inverse
 */

static int check_getri_inverse(PLASMA_desc *descA1, PLASMA_desc *descA2, int *IPIV, double *dparam )
{
    double Rnorm, Anorm, Ainvnorm, result;
    double *work = (double *)malloc(descA1->n*descA1->n*sizeof(double));
    double eps = LAPACKE_dlamch_work('e');
    PLASMA_desc        *descW;

    PLASMA_Desc_Create(&descW, work, PlasmaRealDouble,  descA1->mb, descA1->nb, descA1->bsiz,
                       descA1->m, descA1->n, 0, 0, descA1->m, descA1->n);

    PLASMA_dlaset_Tile( PlasmaUpperLower, (double)0., (double)1., descW);
    PLASMA_dgemm_Tile( PlasmaNoTrans, PlasmaNoTrans,
                       (double)-1., descA2, descA1,
                       (double)1.,  descW);

    Anorm    = PLASMA_dlange_Tile(PlasmaInfNorm, descA1);
    Ainvnorm = PLASMA_dlange_Tile(PlasmaInfNorm, descA2);
    Rnorm    = PLASMA_dlange_Tile(PlasmaInfNorm, descW);

    dparam[IPARAM_ANORM] = Anorm;
    dparam[IPARAM_BNORM] = Ainvnorm;

    result = Rnorm / ( (Anorm*Ainvnorm)*descA1->m*eps ) ;
    dparam[IPARAM_RES] = Rnorm;

    if (  isnan(Ainvnorm) || isinf(Ainvnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        dparam[IPARAM_XNORM] = -1.;
    }
    else{
        dparam[IPARAM_XNORM] = 0.;
    }

    PLASMA_Desc_Destroy(&descW);
    free(work);

    return PLASMA_SUCCESS;
}

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    PLASMA_desc descW;
    int ret = 0;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA,      1, double, PlasmaRealDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA2, check, double, PlasmaRealDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, N, 1 );

    PLASMA_Alloc_Workspace_dgetri_Tile_Async(descA, &descW);
    PLASMA_dplrnt_Tile( descA, 3453 );

    if ( check ) {
        PLASMA_dlacpy_Tile( PlasmaUpperLower, descA, descA2 );
    }

    /* PLASMA DGETRF / DTRTRI / DTRSMRV  */
    {
#if defined(TRACE_BY_SEQUENCE)
        PLASMA_sequence *sequence[4];
        PLASMA_request request[4] = { PLASMA_REQUEST_INITIALIZER,
                                      PLASMA_REQUEST_INITIALIZER,
                                      PLASMA_REQUEST_INITIALIZER,
                                      PLASMA_REQUEST_INITIALIZER };

        PLASMA_Sequence_Create(&sequence[0]);
        PLASMA_Sequence_Create(&sequence[1]);
        PLASMA_Sequence_Create(&sequence[2]);
        PLASMA_Sequence_Create(&sequence[3]);

        if ( ! iparam[IPARAM_ASYNC] ) {

            START_TIMING();
            PLASMA_dgetrf_Tile_Async(descA, piv, sequence[0], &request[0]);
            PLASMA_Sequence_Wait(sequence[0]);

            PLASMA_dtrtri_Tile_Async(PlasmaUpper, PlasmaNonUnit, descA, sequence[1], &request[1]);
            PLASMA_Sequence_Wait(sequence[1]);

            PLASMA_dtrsmrv_Tile_Async(PlasmaRight, PlasmaLower, PlasmaNoTrans, PlasmaUnit,
                                      (double) 1.0, descA, &descW, sequence[2], &request[2]);
            PLASMA_Sequence_Wait(sequence[2]);

            PLASMA_dlaswpc_Tile_Async(descA, 1, descA->m, piv, -1, sequence[3], &request[3]);
            PLASMA_Sequence_Wait(sequence[3]);
            STOP_TIMING();

        } else {

            START_TIMING();
            PLASMA_dgetrf_Tile_Async( descA, piv, sequence[0], &request[0]);
            PLASMA_dtrtri_Tile_Async( PlasmaUpper, PlasmaNonUnit,
                                      descA, sequence[1], &request[1]);
            PLASMA_dtrsmrv_Tile_Async(PlasmaRight, PlasmaLower, PlasmaNoTrans, PlasmaUnit,
                                      (double) 1.0,
                                      descA, &descW, sequence[2], &request[2]);
            PLASMA_dlaswpc_Tile_Async(descA, 1, descA->m, piv, -1,
                                      sequence[3], &request[3]);

            /* Wait for everything */
            PLASMA_Sequence_Wait(sequence[0]);
            PLASMA_Sequence_Wait(sequence[1]);
            PLASMA_Sequence_Wait(sequence[2]);
            PLASMA_Sequence_Wait(sequence[3]);
            STOP_TIMING();

        }

        PLASMA_Sequence_Destroy(sequence[0]);
        PLASMA_Sequence_Destroy(sequence[1]);
        PLASMA_Sequence_Destroy(sequence[2]);
        PLASMA_Sequence_Destroy(sequence[3]);

#else
        if ( ! iparam[IPARAM_ASYNC] ) {

            START_TIMING();
            PLASMA_dgetrf_Tile(descA, piv);
            PLASMA_dtrtri_Tile(PlasmaUpper, PlasmaNonUnit, descA);
            PLASMA_dtrsmrv_Tile(PlasmaRight, PlasmaLower, PlasmaNoTrans, PlasmaUnit,
                                (double) 1.0, descA, &descW);
            PLASMA_dlaswpc_Tile(descA, 1, descA->m, piv, -1);
            STOP_TIMING();

        } else {

            PLASMA_sequence *sequence;
            PLASMA_request request[2] = { PLASMA_REQUEST_INITIALIZER,
                                          PLASMA_REQUEST_INITIALIZER };

            PLASMA_Sequence_Create(&sequence);

            START_TIMING();
            PLASMA_dgetrf_Tile_Async(descA, piv, sequence, &request[0]);
            PLASMA_dgetri_Tile_Async(descA, piv, &descW, sequence, &request[1]);
            PLASMA_Sequence_Wait(sequence);
            STOP_TIMING();

            PLASMA_Sequence_Destroy(sequence);
        }
#endif
    }

    /* Check the solution */
    if ( check )
    {
        ret = check_getri_inverse(descA2, descA, piv, dparam);

        PASTE_CODE_FREE_MATRIX( descA2 );
    }

    PASTE_CODE_FREE_MATRIX( descA );
    free(descW.mat);
    free( piv );

    return ret;
}
