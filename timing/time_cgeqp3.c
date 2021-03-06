/**
 *
 * @generated c Tue Jan  7 11:45:25 2014
 *
 **/
#define _TYPE  PLASMA_Complex32_t
#define _PREC  float
#define _LAMCH LAPACKE_slamch_work

#define _NAME  "PLASMA_cgeqrf"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEQRF(M, N)
#define _FADDS FADDS_GEQRF(M, N)

#include "./timing.c"
#include <cblas.h>

#define COMPLEX

/* computes the residual ||A*P - Q*R||
 * similar to zqpt01.f in testing.
 * A    is m x n
 * AF   is m x n, partial QR factorization of k columns of A.
 * tau  for Householder transformations.
 * jpvt column pivots.
 * work array, dimension (lwork)
 * lwork >= m*n + n.
 */
float qp_residual(
    int m, int n, int k,
    const PLASMA_Complex32_t *A,  int lda,
    const PLASMA_Complex32_t *AF, int ldaf,
    const PLASMA_Complex32_t *tau,
    const int *jpvt,
    PLASMA_Complex32_t *work, int lwork )
{
    PLASMA_Complex32_t mzone = (PLASMA_Complex32_t)-1.0;
    PLASMA_Complex32_t zzero = (PLASMA_Complex32_t) 0.0;
    int j, ldw;
    float rwork[1];

    ldw = m;
    if ( lwork < ldw*n + n ) {
        fprintf( stderr, "Error: lwork %d too small, requires %d\n", lwork, m*n + n );
        return 0.;
    }

    /* copy R to work, with zeros below diagonal */
    LAPACKE_claset( LAPACK_COL_MAJOR, 'L', m-1, n, zzero, zzero, &work[1], ldw );
    LAPACKE_clacpy( LAPACK_COL_MAJOR, 'U', m, n, AF, lda, work, ldw );

    /* form W = QR */
    LAPACKE_cunmqr_work( LAPACK_COL_MAJOR, 'L', 'N', m, n, k,
                         AF, ldaf, tau, work, ldw,
                         &work[n*ldw], lwork-n*ldw );

    //printf( "\n" );
    //printf( "m %d, n %d, k %d\n", m, n, k );
    //CORE_cprint( m, n, A,    lda,  "A"   );
    //CORE_cprint( m, n, AF,   ldaf, "QR" );
    //CORE_cprint( 1, min(m,n), tau, 1, "tau" );
    //CORE_cprint( m, n, work, ldw,  "Q*R" );
    
    /* compare (subtract) j-th column of W to jpvt(j)-th column of A (jpvt is 1-based) */
    for( j = 0; j < n; ++j ) {
        cblas_caxpy( m, CBLAS_SADDR(mzone), &A[(jpvt[j]-1)*lda], 1, &work[j*ldw], 1 );
    }
    
    //CORE_cprint( m, n, work, ldw,  "Q*R - AP" );
    
    return LAPACKE_clange_work( LAPACK_COL_MAJOR, 'O', m, n, work, ldw, rwork );
}


static int
RunTest(int *iparam, float *dparam, real_Double_t *t_)
{
    int info, lwork, lwork2, i;
    int *jpvt;
    PLASMA_Complex32_t *tau, *work, *work2;
    float *rwork;

    PASTE_CODE_IPARAM_LOCALS( iparam );

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, PLASMA_Complex32_t, LDA, N );

    /* LAPACK needs larger size in real;
     * PLASMA currently uses work and rwork for both real and complex. */
    #ifdef COMPLEX
    lwork = (N+1)*NB;
    #else
    lwork = (N+1)*NB + 2*N;
    #endif
    work  = (PLASMA_Complex32_t*) malloc( lwork * sizeof(PLASMA_Complex32_t) );
    rwork = (float*)             malloc( 2*2*N * sizeof(float)             );
    jpvt  = (int*)                malloc( N     * sizeof(int)                );
    tau   = (PLASMA_Complex32_t*) malloc( 2*N   * sizeof(PLASMA_Complex32_t) );

    if ( jpvt == NULL || tau == NULL || work == NULL || rwork == NULL ) {
        fprintf( stderr, "malloc failed\n" );
        return -1;
    }

    /* zero out pivots (required by LAPACK) */
    for( i = 0; i < N; ++i ) {
        jpvt[i] = 0;
    }

    /* Initialize Data */
    PLASMA_cplrnt(M, N, A, LDA, 123456);

    /* Save A in lapack layout for check */
    PASTE_CODE_ALLOCATE_COPY( Acpy, check, PLASMA_Complex32_t, A, LDA, N );

    START_TIMING();
    info = PLASMA_cgeqp3( M, N, A, LDA, jpvt, tau, work, rwork );
    STOP_TIMING();

    /* Check the solution */
    if ( info != 0 ) {
        printf( "\nPLASMA_cgeqp3 returned error %d.\n", info );
    }
    else if ( check ) {
        lwork2 = (M*N + N);
        work2  = (PLASMA_Complex32_t*) malloc( lwork2 * sizeof(PLASMA_Complex32_t) );
        if ( work2 == NULL ) {
            fprintf( stderr, "test malloc failed\n" );
            return -1;
        }

        dparam[IPARAM_ANORM] = LAPACKE_clange_work( LAPACK_COL_MAJOR, 'F', N, N, Acpy, LDA, rwork );
        dparam[IPARAM_XNORM] = 1.;
        dparam[IPARAM_BNORM] = 0.;
        dparam[IPARAM_RES]   = qp_residual( M, N, min(M,N), Acpy, LDA,
                                            A, LDA, tau, jpvt, work2, lwork2 );

        // /* compute result with LAPACK */
        // int *jpvt2;
        // PLASMA_Complex32_t *tau2;
        // jpvt2  = (int*)                malloc( N * sizeof(int)    );
        // tau2   = (PLASMA_Complex32_t*) malloc( N * sizeof(PLASMA_Complex32_t) );
        // if ( jpvt2 == NULL || tau2 == NULL ) {
        //     fprintf( stderr, "test malloc failed\n" );
        //     return -1;
        // }
        // 
        // /* zero out pivots (required by LAPACK) */
        // for( i = 0; i < N; ++i ) {
        //     jpvt2[i] = 0;
        // }
        // 
        // float time = cWtime();
        // #ifdef COMPLEX
        // info = LAPACKE_cgeqp3_work( LAPACK_COL_MAJOR, N, N, Acpy, LDA, jpvt2, tau2, work, lwork, rwork );
        // #else
        // info = LAPACKE_cgeqp3_work( LAPACK_COL_MAJOR, N, N, Acpy, LDA, jpvt2, tau2, work, lwork );
        // #endif
        // time = cWtime() - time;
        // if ( info != 0 ) {
        //     printf( "qp3 returned error %d\n", info );
        // }
        // /* printf( "   %7.3f", time ); */
        //         
        // CORE_cprint( M, N, Acpy, LDA, "QR_L" );
        // CORE_cprint( 1, min(M,N), tau2, 1, "tau_L" );
        
        
        // cblas_caxpy( LDA*N, CBLAS_SADDR(mzone), A, 1, Acpy, 1 );
        // dparam[IPARAM_RES] = LAPACKE_clange_work( LAPACK_COL_MAJOR, 'F', N, N, Acpy, LDA, rwork );
        // 
        // cblas_caxpy( N, CBLAS_SADDR(mzone), tau, 1, tau2, 1 );
        // tnorm = cblas_scnrm2( N, tau2, 1 );
        /* printf( "|t|=%.2e\n", tnorm ); */
        
        // for( i = 0; i < N; ++i ) {
        //     if ( jpvt[i]+1 != jpvt2[i] ) {
        //         printf( "pivot mis-match jpvt[%2d]+1=%2d, jpvt2[%2d]=%2d\n", i, jpvt[i]+1, i, jpvt2[i] );
        //     }
        // }
        
        // free( jpvt2 );
        // free( tau2  );

        free( Acpy  );
        free( work2 );
    }

    /* Free data */
    free( A     );
    free( jpvt  );
    free( tau   );
    free( work  );
    free( rwork );

    return 0;
}
