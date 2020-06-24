/**
 * @file dsyev.c
 *
 *  PLASMA computational routines
 *  Release Date: November, 15th 2009
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Azzam Haidar
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @generated d Tue Jan  7 11:45:09 2014
 *
 **/
#include <lapacke.h>
#include "common.h"
/***************************************************************************//**
 *
 * @ingroup double
 *
 *  PLASMA_dsyev - Computes all eigenvalues and, optionally,
 *  eigenvectors of a complex Hermitian matrix A. The matrix A is
 *  preliminary reduced to tridiagonal form using a two-stage
 *  approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If uplo = PlasmaLower, the leading N-by-N lower triangular
 *          part of A contains the lower triangular part of the matrix
 *          A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, the lower triangle (if uplo = PlasmaLower) or the
 *          upper triangle (if uplo = PlasmaUpper) of A, including the
 *          diagonal, is destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_dsyev
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsyev_Tile
 * @sa PLASMA_dsyev_Tile_Async
 * @sa PLASMA_cheev
 * @sa PLASMA_dsyev
 * @sa PLASMA_ssyev
 *
 ******************************************************************************/
int PLASMA_dsyev(PLASMA_enum jobz, PLASMA_enum uplo, int N,
                 double *A, int LDA,
                 double *W,
                 PLASMA_desc *descT,
                 double *Q, int LDQ)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA_dsyev", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_dsyev", "illegal value of jobz");
        return -1;
    }
    if (uplo != PlasmaLower && uplo != PlasmaUpper) {
        plasma_error("PLASMA_dsyev", "illegal value of uplo");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_dsyev", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_dsyev", "illegal value of LDA");
        return -5;
    }
    if (LDQ < max(1, N)) {
        plasma_error("PLASMA_dsyev", "illegal value of LDQ");
        return -9;
    }

    /* Quick return */
    if (N == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_DSYEV, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsyev", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_diplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_dsyev_Tile_Async(jobz, uplo, &descA, W, descT, Q, LDQ, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_dooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_diptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}
/***************************************************************************//**
 *
 * @ingroup double_Tile
 *
 *  PLASMA_dsyev_Tile - Computes all eigenvalues and, optionally, eigenvectors of a
 *  complex Hermitian matrix A using a two-stage approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          A contains the lower triangular part of the matrix A, and
 *          the strictly upper triangular part of A is not referenced.
 *          On exit, if jobz = PlasmaVec, then if return value = 0, A
 *          contains the orthonormal eigenvectors of the matrix A.
 *          If jobz = PlasmaNoVec, then on exit the lower triangle (if
 *          uplo = PlasmaLower) or the upper triangle (if uplo =
 *          PlasmaUpper) of A, including the diagonal, is destroyed.*
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in,out] T
 *          On entry, descriptor as return by
 *          PLASMA_Alloc_Workspace_dsyev
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 * @param[in] LDQ
 *          The leading dimention of the eigenvectors matrix Q. LDQ >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if INFO = i, the algorithm failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsyev_Tile
 * @sa PLASMA_dsyev_Tile_Async
 * @sa PLASMA_cheev
 * @sa PLASMA_dsyev
 * @sa PLASMA_ssyev
 *
 ******************************************************************************/
int PLASMA_dsyev_Tile(PLASMA_enum jobz, PLASMA_enum uplo,
                      PLASMA_desc *A, double *W,
                      PLASMA_desc *T, double *Q, int LDQ)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsyev_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_dsyev_Tile_Async(jobz, uplo, A, W, T, Q, LDQ, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup double_Tile_Async
 *
 *  PLASMA_dsyev_Tile_Async - Computes all eigenvalues and,
 *  optionally, eigenvectors of a complex Hermitian matrix A using a
 *  two-stage approach:
 *  First stage: reduction to band tridiagonal form;
 *  Second stage: reduction from band to tridiagonal form.
 *
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of
 *          A contains the lower triangular part of the matrix A, and
 *          the strictly upper triangular part of A is not referenced.
 *          On exit, if jobz = PlasmaVec, then if return value = 0, A
 *          contains the orthonormal eigenvectors of the matrix A.
 *          If jobz = PlasmaNoVec, then on exit the lower triangle (if
 *          uplo = PlasmaLower) or the upper triangle (if uplo =
 *          PlasmaUpper) of A, including the diagonal, is destroyed.*
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in,out] T
 *          On entry, descriptor as return by
 *          PLASMA_Alloc_Workspace_dsyev
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 * @param[in] LDQ
 *          The leading dimention of the eigenvectors matrix Q. LDQ >= max(1,N).
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa PLASMA_dsyev
 * @sa PLASMA_dsyev_Tile
 * @sa PLASMA_cheev_Tile_Async
 * @sa PLASMA_dsyev_Tile_Async
 * @sa PLASMA_ssyev_Tile_Async
 *
 ******************************************************************************/
int PLASMA_dsyev_Tile_Async(PLASMA_enum jobz, PLASMA_enum uplo,
                            PLASMA_desc *A,
                            double *W,
                            PLASMA_desc *T,
                            double *Q, int LDQ,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA;
    PLASMA_desc descT;
    PLASMA_desc descQ;
    double *AB;
    double *E;
    int N;
    int NB;
    int LDAB;
    int i;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_dsyev_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_dsyev_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_dsyev_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsyev_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsyev_Tile_Async", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }

    /* Check input arguments */
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_dsyev_Tile_Async", "illegal value of jobz");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaLower && uplo != PlasmaUpper) {
        plasma_error("PLASMA_dsyev_Tile_Async", "illegal value of uplo");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descA.m != descA.n) {
        plasma_error("PLASMA_dsyev_Tile_Async", "matrix need to be square");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_dsyev_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    N  = descA.m;
    NB = min(descA.mb,descA.m);
    LDAB = 2*NB+1;

    /* Allocate workspace for band storage of the band matrix A and for the off diagonal after tridiagonalisation */
    AB = (double *)plasma_shared_alloc(plasma, LDAB*N, PlasmaRealDouble);
    memset( AB, 0, LDAB * N * sizeof(double) );
    if (AB == NULL) {
        plasma_error("PLASMA_dsyev_Tile_Async", "plasma_shared_alloc(AB) failed");
        plasma_shared_free(plasma, AB);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    E = (double *)plasma_shared_alloc(plasma, N, PlasmaRealDouble);
    if (E == NULL) {
        plasma_error("PLASMA_dsyev_Tile_Async", "plasma_shared_alloc(E) failed");
        plasma_shared_free(plasma, E);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    #if defined(ENABLE_TIMER)
    PLASMA_Double_t timelpk=0.0,timeaplQ2=0.0, timeT=0.0;
    PLASMA_Double_t timeB=0.0,timeblg=0.0,timeaplQ1=0.0,timeconv1=0.0,timeconv2=0.0,timeall=0.0;
    timeall = PLASMA_Wtime();
    timeB   = PLASMA_Wtime();
    #endif
    /* Reduction to tridiagonal form
     * with a two-stage approach.
     */
    /*=======================================
     *  calling Reduction to BAND
     *  then convert matrix to band form
     *=======================================*/
    plasma_dynamic_call_5(plasma_pdsytrd_he2hb,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_call_6( plasma_pdsbcpy_t2bl,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        double*, AB,
        int, LDAB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_sync();

    status = sequence->status;
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_dsyev","pdsytrd_he2hb+pzcopy");
        return status;
    }
    /*=======================================
     *  END of calling Reduction to BAND
     *=======================================*/
    #if defined(ENABLE_TIMER)
    timeB   = PLASMA_Wtime()-timeB;
    printf("\n  Finish Band red    timing= %lf \n",timeB);
    #endif
    /*=======================================
     *  calling bulge chasing
     *=======================================*/
    double *TAU2 = NULL;
    double *V2 = NULL;
    double *T2 = NULL;
    int Vblksiz, blkcnt, LDT, LDV;
    int WANTZ   = 0;
    int blguplo = PlasmaLower;
    /* int NE      = N; // for later use when a portion of the eigenvectors are requested*/
    if( jobz == PlasmaNoVec )
        WANTZ=0;
    else
        WANTZ=1;

    /* Vblksiz correspond to the blocking used when applying V2 to the matrix Q
     * it is similar to IB in LAPACK DORMQR.
     * blkcnt is the number of losange or tile of Vs */
    /* Note that in case PlamaVec requested, the V2 and T2 are stored by the
     * bulgechasing function in a special format:
     * for V2s: it store the V2(LDV,Vblksiz) of each losange in a tile storage meaning
     * that V2_1 is stored then V2_2,..., V2_blkcnt.
     * blkcnt is the number of losange.
     * */
    Vblksiz = min(NB,48);
    LDT     = Vblksiz;
    if( jobz == PlasmaVec ) {
        findVTsiz(N, NB, Vblksiz, &blkcnt, &LDV);
        TAU2   = (double *) plasma_shared_alloc(plasma,     blkcnt*Vblksiz, PlasmaRealDouble);
        V2     = (double *) plasma_shared_alloc(plasma, LDV*blkcnt*Vblksiz, PlasmaRealDouble);
        T2     = (double *) plasma_shared_alloc(plasma, LDT*blkcnt*Vblksiz, PlasmaRealDouble);
        if ( (TAU2 == NULL) || (V2 == NULL) || (T2 == NULL) ) {
            plasma_error("PLASMA_dsyev", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAU2);
            plasma_shared_free(plasma, V2);
            plasma_shared_free(plasma, T2);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAU2, 0,     blkcnt*Vblksiz*sizeof(double));
        memset(V2,   0, LDV*blkcnt*Vblksiz*sizeof(double));
        memset(T2,   0, LDT*blkcnt*Vblksiz*sizeof(double));
    }
    else {
        TAU2   = (double *) plasma_shared_alloc(plasma, 2*N, PlasmaRealDouble);
        V2     = (double *) plasma_shared_alloc(plasma, 2*N, PlasmaRealDouble);
        if ( (TAU2 == NULL) || (V2 == NULL) ) {
            plasma_error("PLASMA_dsyev", "plasma_shared_alloc() failed");
            plasma_shared_free(plasma, TAU2);
            plasma_shared_free(plasma, V2);
            return PLASMA_ERR_OUT_OF_RESOURCES;
        }
        memset(TAU2, 0, 2*N*sizeof(double));
        memset(V2,   0, 2*N*sizeof(double));
    }
    #if defined(ENABLE_TIMER)
    timeblg  = PLASMA_Wtime();
    #endif
    plasma_parallel_call_13(plasma_pdsytrd_hb2st_v1,
        PLASMA_enum,         blguplo,
        int,                 N,
        int,                 NB,
        int,                 Vblksiz,
        double*, AB,
        int,                 LDAB,
        double*, V2,
        double*, TAU2,
        double*,             W,
        double*,             E,
        int,                 WANTZ,
        PLASMA_sequence*,    sequence,
        PLASMA_request*,     request);
    /* WARNING: If plasma_pdsytrd_hb2st is implemented through a dynamic call, don't
     * forget to synchronize */
    plasma_dynamic_sync();
    /*=======================================
     *  END of calling bulge chasing
     *=======================================*/
    #if defined(ENABLE_TIMER)
    timeblg   = PLASMA_Wtime()-timeblg;
    printf("  Finish Bulge       timing= %lf \n" ,timeblg);
    #endif
    if (jobz == PlasmaVec){
        /*=======================================
         *  generate Q2 from the bulge
         *=======================================*/
        /* Initialize Q to Identity */
        memset(Q,   0, N*LDQ*sizeof(double));
        for(i=0; i<N; i++){
            Q[i+i*LDQ] = 1.0;
        }
        /* compute T2 */
        #if defined(ENABLE_TIMER)
        timeT  = PLASMA_Wtime();
        #endif
        plasma_static_call_8(plasma_pdlarft_blgtrd,
            int,                 N,
            int,                 NB,
            int,                 Vblksiz,
            double*, V2,
            double*, T2,
            double*, TAU2,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);
        #if defined(ENABLE_TIMER)
        plasma_dynamic_sync();
        timeT  = PLASMA_Wtime()-timeT;
        printf("  Finish compute T2  timing= %lf \n" ,timeT);
        timeaplQ2  = PLASMA_Wtime();
        #endif
        /* apply Q2 from Left */
        plasma_static_call_14(plasma_pdormqr_blgtrd,
            PLASMA_enum,         PlasmaLeft,
            PLASMA_enum,         PlasmaNoTrans,
            int,                 N,
            int,                 NB,
            int,                 N,
            int,                 Vblksiz,
            int,                 WANTZ,
            double*, V2,
            double*, T2,
            double*, TAU2,
            double*, Q,
            int,                 LDQ,
            PLASMA_sequence*,    sequence,
            PLASMA_request*,     request);
        #if defined(ENABLE_TIMER)
        plasma_dynamic_sync();
        timeaplQ2  = PLASMA_Wtime()-timeaplQ2;
        printf("  Finish compute Q2  timing= %lf \n" ,timeaplQ2);
        #endif
        /*=======================================
         *  apply Q1 from the reduction to band
         *=======================================*/
        /* CASE NB>N, Q1 doesn't need to be applied, only bulge chasing has been done */
        if( NB < N ){
            #if defined(ENABLE_TIMER)
            plasma_dynamic_sync();
            timeconv1   = PLASMA_Wtime();
            #endif
            if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
                plasma_dooplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, N, N, sequence, request, plasma_desc_mat_free(&(descQ)) );
            } else {
                plasma_diplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, N, N, sequence, request);
            }
            #if defined(ENABLE_TIMER)
            timeconv1    = PLASMA_Wtime()-timeconv1;
            timeaplQ1   = PLASMA_Wtime();
            #endif
            /* Accumulate the transformations from the first stage */
            if(uplo==PlasmaLower){
                plasma_dynamic_call_7(plasma_pdormqr,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaNoTrans,
                    PLASMA_desc, plasma_desc_submatrix(descA, descA.mb, 0, descA.m-descA.mb, descA.n-descA.nb),
                    PLASMA_desc, plasma_desc_submatrix(descQ, descQ.mb, 0, descQ.m-descQ.mb, descQ.n),
                    PLASMA_desc, plasma_desc_submatrix(descT, descT.mb, 0, descT.m-descT.mb, descT.n-descT.nb),
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);

            }
            else {
                plasma_dynamic_call_7(plasma_pdormlq,
                    PLASMA_enum, PlasmaLeft,
                    PLASMA_enum, PlasmaTrans,
                    PLASMA_desc, plasma_desc_submatrix(descA, 0, descA.nb, descA.m-descA.mb, descA.n-descA.nb),
                    PLASMA_desc, plasma_desc_submatrix(descQ, descQ.mb, 0, descQ.m-descQ.mb, descQ.n),
                    PLASMA_desc, plasma_desc_submatrix(descT, 0, descT.nb, descT.m-descT.mb, descT.n-descT.nb),
                    PLASMA_sequence*, sequence,
                    PLASMA_request*, request);
            }
            #if defined(ENABLE_TIMER)
            plasma_dynamic_sync();
            timeaplQ1   = PLASMA_Wtime()-timeaplQ1;
            printf("  Finish compute Q1  timing= %lf \n" ,timeaplQ1);
            timeconv2    = PLASMA_Wtime();
            #endif
            if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
                plasma_dooptile2lap( descQ, Q, NB, NB, LDQ, N, sequence, request );
                plasma_dynamic_sync();
                plasma_desc_mat_free(&descQ);
            } else {
                plasma_diptile2lap( descQ, Q, NB, NB, LDQ, N, sequence, request );
                plasma_dynamic_sync();
            }
            #if defined(ENABLE_TIMER)
            plasma_dynamic_sync();
            timeconv2    = PLASMA_Wtime()-timeconv2;
            printf("  Finish convert     timing= %lf \n" ,timeconv1+timeconv2);
            #endif
        } /* END of ( NB < N ) */
    }
    /*=======================================
     *  END of calling computing Q
     *=======================================*/
    #if defined(ENABLE_TIMER)
    plasma_dynamic_sync();
    timelpk   = PLASMA_Wtime();
    #endif
    /*=======================================
     *  calling eigensolver
     *=======================================*/
    plasma_setlapack_multithreads(plasma->world_size);
    /* call eigensolver using lapack routine for our resulting tridiag [D E] */
    status = LAPACKE_dsteqr( LAPACK_COL_MAJOR, lapack_const(jobz),
                    N, W, E, Q, LDQ );
    if(status != 0){
        plasma_error("PLASMA_dstedc","ZSTEQR");
        /*return status;*/
    }
    sequence->status = status;
    plasma_setlapack_sequential(plasma);
    /*=======================================
     *  END of calling eigensolver
     *=======================================*/
    #if defined(ENABLE_TIMER)
    timelpk = PLASMA_Wtime()-timelpk;
    printf("  Finish Eigensolver timing= %lf  WANTZ %d threads %d\n" ,timelpk, WANTZ,  plasma->world_size);
    timeall = PLASMA_Wtime()-timeall;
    printf("  Finish full eigenproblem threads %d  N %d  timeall= %lf \n", plasma->world_size, N, timeall);
    #endif

    if( jobz == PlasmaVec )
        plasma_shared_free(plasma, T2);
    plasma_shared_free(plasma, V2);
    plasma_shared_free(plasma, TAU2);
    plasma_shared_free(plasma, E);
    plasma_shared_free(plasma, AB);
    return PLASMA_SUCCESS;
}
